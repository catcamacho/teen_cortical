#! /bin/csh

set raw_dir = /share/iang/active/ELS/ELS_FreeSurfer/ELS_FS_subjDir
setenv SUBJECTS_DIR /share/iang/active/ELS/ELS_FreeSurfer/Analysis/proc/template_fs6

foreach sub(145-T1 156-Tmid 137x-Tmid 183-T1 307-TK1 031-T1 196-T1 202-T1 164-Tmid 005-T1 125-Tmid 146-Tmid 002-T1 034-T1 159-T2 302-TK1 026-T1 120-T2 023-T1 110x-T2 069-T2 086x-T1 103-T1 106-T2 055-T2 077-T1 149-T1 024-Tmid 038-T2 006-TK3)
    set raw = ${raw_dir}/${sub}/mri/orig/001.mgz
    recon-all -s ${sub} -i $raw -all -gcut -openmp 8
    
    if(-e $SUBJECTS_DIR/${sub}/scripts/recon-all.done) then
        echo "--------------" ${sub} " finished running." >> $SUBJECTS_DIR/log.txt
    else
        echo "WARNING: " ${sub} " did not complete!!" >> $SUBJECTS_DIR/log.txt
    endif
end

