
# coding: utf-8

# # ELS Cortical Thickness Processing Workflows
# This notebook is split into several parts:
# * Study-wide variables
# * Template creation workflow
# * Template subject tissue segmentation workflow
# * Tissue priors workflow
# * Subject cortical thickness workflow
# 
# Template subject segmentations are edited before being used to create the tissue priors.

# In[ ]:


# Import modules
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink, FreeSurferSource
from nipype.interfaces.fsl.preprocess import FAST
from nipype.interfaces.fsl.utils import Reorient2Std
from nipype.interfaces.freesurfer import FSCommand, MRIConvert, Binarize
from os import listdir

# Set up study specific variables
project_home = '/share/iang/active/ELS/ELS_FreeSurfer/Analysis'
#project_home = '/Users/myelin/Dropbox/data/fs_practice'
#project_home = '/Users/lucindasisk/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice'
#project_home = '/Users/catcamacho/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice/Analysis'

#fs_subjdir = '/share/iang/active/ELS/ELS_FreeSurfer/ELS_FS_subjDir'
#fs_subjdir = '/Users/myelin/Dropbox/data/fs_practice/ELS_FS_subjDir'
#fs_subjdir = '/Users/catcamacho/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice/ELS_FS_subjDir'
#fs_subjdir = '/Users/lucindasisk/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice/ELS_FS_subjDir'
fs_subjdir = '/share/iang/active/ELS/ELS_FreeSurfer/Analysis/proc/template_fs6' 

workflow_dir = project_home + '/workflows'
subj_proc = project_home + '/proc/subject'
group_proc = project_home + '/proc/group'
template_proc = project_home + '/proc/template'
#subject_info = project_home + '/misc/subjects.csv' 
#template_sub = ['011-T1']
template_sub=listdir(fs_subjdir)

#set default FreeSurfer subjects dir
FSCommand.set_default_subjects_dir(fs_subjdir)


# In[ ]:


######### File handling #########

#Pass in list to freesurfer source node (subs) 
fs_source = MapNode(FreeSurferSource(subjects_dir = fs_subjdir), 
                    name = 'fs_source', iterfield = ['subject_id'])
fs_source.inputs.subject_id = template_sub

#set up datasink
datasink = Node(DataSink(base_directory = template_proc),
                name = 'datasink')


# ## Template Creation Workflow
# Below are the cells associated with template creation:
# * Unique functions
#     - make3DTemplate wraps the ANTs antsMultivariateTemplateConstruction2 script
# * Template creation Nodes
# * Template creation workflow steps
#     - Convert template subjects' FreeSurfer T1 images to nifti 
#     - Reorient T1s to standard
#     - Pass the T1s to the ANTs template creation script

# In[ ]:


######### Template creation functions #########
def make3DTemplate(subject_T1s, num_proc, output_prefix):
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    from os.path import abspath, split
    from os import getcwd
    from shutil import copyfile
    from glob import glob
    from subprocess import call

    curr_dir = getcwd()

    #copy T1s into current directory
    for T in range(0,len(subject_T1s)):
        [dirname,filename] = split(subject_T1s[T])
        copyfile(subject_T1s[T],curr_dir + '/S' + str(T)+'_'+filename)

    # -c flag is control for local computing (2= use localhost; required for -j flag)
    # -j flag is for number of processors allowed
    call(['antsMultivariateTemplateConstruction2.sh', '–d','3','–o', output_prefix,'–r','1','–c','2','–j', str(num_proc), '*.nii.gz'])
    
    sample_template = abspath(output_prefix + 'template0.nii.gz')
    
    return(sample_template)


# In[ ]:


######### Template creation nodes #########

#convert freesurfer brainmask files to .nii
convertT1 = MapNode(MRIConvert(out_file='brainmask.nii.gz',
                               out_type='niigz'), 
                    name='convertT1', 
                    iterfield = ['in_file'])

#reorient files to standard space
reorientT1 = MapNode(Reorient2Std(),
                     name = 'reorientT1',
                     iterfield = ['in_file'])

#pass files into template function (normalized, pre-skull-stripping)
makeTemplate = Node(Function(input_names=['subject_T1s','num_proc','output_prefix'],
                             output_names=['sample_template'],
                             function=make3DTemplate),
                    name='makeTemplate')
makeTemplate.inputs.num_proc=8 # feel free to change to suit what's free on SNI-VCS
makeTemplate.inputs.output_prefix='ELS_CT_'


# In[ ]:


######### Template creation workflow #########
template_flow = Workflow(name = "template_flow")
template_flow.connect([(fs_source, convertT1, [('T1','in_file')]),
                       (convertT1, reorientT1, [('out_file', 'in_file')]),
                       #(reorientT1, makeTemplate, [('out_file', 'subject_T1s')]),
                       (reorientT1, datasink, ['out_file','sub_T1'])
                       #(makeTemplate, datasink, [('sample_template', 'sample_template')])
                      ])

template_flow.base_dir = workflow_dir
template_flow.write_graph(graph2use = 'flat')
template_flow.run()


# ## Template Subject Tissue Segmentation Workflow
# These cells are associated with the template subjects tissue segmentation workflow. The point is ultimately create 5 tissue class per subject: cerebrospinal fluid, cortical gray matter, subcortical gray matter, white matter, and whole brain.
# * Custom functions
# * Nodes
# * Workflow
#     - FSL's FAST is used to segment CSF and WM
#     - FreeSurfer's automated segmentation is used to segment subcortical and cortical GM
#     - FreeSurfer's brainmask is binarized for whole brain extraction
# 
# The resulting masks are manually edited for accuracy.

# In[ ]:


######### Tissue creation nodes: segment template subjects to 5 tissue classes #########
### 1)CSF, 2)cortical gray matter, 3)white matter, 4)subcortical gray matter, and 5)whole brain

def relabel_fast(fast_tissue_list):
    from nipype import config, logging
    from os.path import split
    from os import rename
    config.enable_debug_mode()
    logging.update_logging(config)
    tissue_list = sorted(fast_tissue_list)
    csf = tissue_list[0]
    wm = tissue_list[2]
    [wd, csf_file] = split(csf)
    [wd, wm_file] = split(wm)
    rename(csf, wd + 'csf.nii.gz')
    rename(wm, wd + 'wm.nii.gz')
    wm_csf = [wd + 'csf.nii.gz', wd + 'wm.nii.gz']
    return(wm_csf)

# Create subcortical and cortical gray matter masks <-- custom function. inputs: fs aseg + tissue class files
def aseg_to_tissuemaps(aseg):
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    from nibabel import load, save, Nifti1Image
    from numpy import zeros_like
    from os.path import abspath
    aseg_nifti = load(aseg)
    aseg_data = aseg_nifti.get_data()
    cortical_labels = [3, 42]
    subcortical_labels =[8, 10, 11, 12, 13, 17, 18, 26, 47, 49, 50, 51, 52, 53, 54, 58]

    #creating array of zeroes that replaces 0's with 1's when matches values of subcortical_labels
    cortical_data = zeros_like(aseg_data)
    for x in cortical_labels:
        cortical_data[aseg_data == x] = 1
    cortical_nifti = Nifti1Image(cortical_data, aseg_nifti.affine)
    
    subcort_data = zeros_like(aseg_data) 
    for x in subcortical_labels:
        subcort_data[aseg_data == x] = 1
    subcort_nifti = Nifti1Image(subcort_data, aseg_nifti.affine)
    
    save(subcort_nifti, "subcortical_gm.nii.gz")
    save(cortical_nifti, "cortical_gm.nii.gz")
    subcort_file = abspath("subcortical_gm.nii.gz")
    cortical_file = abspath("cortical_gm.nii.gz")
    gm_list = [subcort_file, cortical_file]
    return(gm_list)


# In[ ]:


#convert freesurfer brainmask files to .nii
convert_to_nii = MapNode(MRIConvert(out_file='brainmask.nii.gz',
                                    out_type='niigz'), 
                         name='convert_to_nii', 
                         iterfield = ['in_file'])

#reorient brainmask file to standard
reorient_to_std = MapNode(Reorient2Std(),
                         name = 'reorient_to_std',
                         iterfield = ['in_file'])

# Reorient aseg to standard
reorient_aseg = MapNode(Reorient2Std(),
                         name = 'reorient_aseg',
                         iterfield = ['in_file'])

#T1 gets run through segmentation (2) ---> results in segmentation into 3 tissue classes (wm, gm, csf)
segment = MapNode(FAST(number_classes = 3, 
                       segments=True, 
                       no_bias=True), 
                  name = 'segment', 
                  iterfield = ['in_files'])

# Convert freesurfer aseg to nii
convert_aseg = MapNode(MRIConvert(out_file='aseg.nii.gz',
                                    out_type='niigz'), 
                         name='convert_aseg', 
                         iterfield = ['in_file'])

# Split aseg into to types of gray matter
aseg_to_gm = MapNode(Function(input_names=['aseg'],
                              output_names=['gm_list'],
                              function=aseg_to_tissuemaps),
                     name='aseg_to_gm', 
                     iterfield=['aseg'])

# Relabel the FAST segmentation 
relabel_fast_seg = MapNode(Function(input_names=['fast_tissue_list'],
                                    output_names=['wm_csf'],
                                    function=relabel_fast),
                           name='relabel_fast_seg',
                           iterfield=['fast_tissue_list'])

# make brainmask by binarizing the brainmask
binarize_brain = MapNode(Binarize(min=1, 
                                  dilate=1, 
                                  erode=1), 
                         name='binarize_brain', 
                         iterfield=['in_file'])


# In[ ]:


######### Tissue segmentation workflow #########
segment_flow = Workflow(name = "segment_flow")
segment_flow.connect([(fs_source, convert_to_nii, [('brainmask','in_file')]),
                      (convert_to_nii, reorient_to_std, [('out_file', 'in_file')]),
                      (reorient_to_std, segment, [('out_file', 'in_files')]),
                      (segment, relabel_fast_seg, [('tissue_class_files', 'fast_tissue_list')]),
                      (fs_source, convert_aseg, [('aseg','in_file')]),
                      (convert_aseg, reorient_aseg, [('out_file', 'in_file')]),
                      (reorient_aseg, aseg_to_gm, [('out_file', 'aseg')]),
                      (reorient_to_std, binarize_brain, [('out_file','in_file')]),
                      (binarize_brain, datasink, [('binary_file','brain_seg')]),
                      (aseg_to_gm, datasink, [('gm_list', 'gm_files')]),
                      (relabel_fast_seg, datasink, [('wm_csf', 'wm_csf')])
                     ])

segment_flow.base_dir = workflow_dir
segment_flow.write_graph(graph2use = 'flat')
segment_flow.run('MultiProc', plugin_args={'n_procs': 10})