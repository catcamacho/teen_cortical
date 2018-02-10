
# coding: utf-8

# In[1]:


# Import modules
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink, FreeSurferSource
from nipype.interfaces.fsl.preprocess import FAST
from nipype.interfaces.fsl.utils import Reorient2Std
from nipype.interfaces.freesurfer import FSCommand, MRIConvert
from os import listdir

# Set up study specific variables
#project_home = '/Volumes/iang/active/ELS/ELS_FreeSurfer/Analysis'
#project_home = '/Users/myelin/Dropbox/data/fs_practice'
#project_home = '/Users/catcamacho/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice'
project_home = '/share/iang/active/ELS/ELS_FreeSurfer/Analysis'

#fs_subjdir = '/Volumes/iang/active/ELS/ELS_FreeSurfer/ELS_FS_subjDir'
#fs_subjdir = '/Users/myelin/Dropbox/data/fs_practice/ELS_FS_subjDir'
#fs_subjdir = '/Users/catcamacho/Dropbox/Projects/AC_ELS_Cortex/proc/fs_practice/ELS_FS_subjDir'
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


# In[2]:


######### File handling #########

#Pass in list to freesurfer source node (subs) 
fs_source = MapNode(FreeSurferSource(subjects_dir = fs_subjdir), 
                    name = 'fs_source', iterfield = ['subject_id'])
fs_source.inputs.subject_id = template_sub

#set up datasink
datasink = Node(DataSink(base_directory = template_proc),
                name = 'datasink')


# In[3]:


######### Template creation functions #########
def make3DTemplate(subject_T1s, num_proc, output_prefix):
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    from os.path import abspath, split
    from os import getcwd
    from shutil import copyfile
    from glob import glob

    curr_dir = getcwd()

    #copy T1s into current directory
    for T1 in subject_T1s:
        copyfile(T1,curr_dir)

    # determine the common suffix across all the files in the subject_T1s list
    T1_suffix = subject_T1s[0]
    for a in range(0,len(subject_T1s)):
        while T1_suffix not in subject_T1s[a]:
            T1_suffix = T1_suffix[1:]
    
    # determine the common prefix across all the files in the subject_T1s list
    (folder, T1_prefix) = split(subject_T1s[0])
    for a in range(0,len(subject_T1s)):
        while T1_prefix not in subject_T1s[a]:
            T1_prefix = T1_prefix[:-1]

    # -c flag is control for local computing (2= use localhost; required for -j flag)
    # -j flag is for number of processors allowed
    call(['antsMultivariateTemplateConstruction2.sh –d 3 –o %s –r 1 –c 2 –j %d %s*%s' % (output_prefix, num_proc, T1_prefix, T1_suffix)])
    
    sample_template = abspath(output_prefix + 'template0.nii.gz')
    
    return(sample_template)


# In[4]:


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
makeTemplate.inputs.num_proc=16 # feel free to change to suit what's free on SNI=VCS
makeTemplate.inputs.output_prefix='ELS_CT_'


# In[5]:


######### Template creation workflow #########
template_flow = Workflow(name = 'template_flow')
template_flow.connect([(fs_source, convertT1, [('T1','in_file')]),
                       (convertT1, reorientT1, [('out_file', 'in_file')]),
                       (reorientT1, makeTemplate, [('out_file', 'subject_T1s')]),
                       (makeTemplate, datasink, [('sample_template', 'sample_template')])
                      ])

template_flow.base_dir = workflow_dir
template_flow.write_graph(graph2use = 'flat')
template_flow.run()

