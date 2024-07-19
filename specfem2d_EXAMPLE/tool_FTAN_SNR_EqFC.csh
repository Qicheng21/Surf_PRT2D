

#foreach event (20200331235231)
#foreach event (20200318130931 20200331235231)
#cd $event
#ls $event"".*EHZ*sac | cut -d. -f2 > stalst   
#lf_touch_sac_rotation_GCP stalst $event
ls *EHZ.sac *EHR.sac > sac_R.lst
#lf_touch_sac_v1 sac_R.lst

#ls *EHZ.sac *EHR.sac | awk '{print "0 1 4.5 0.1 30 20 1 0.5 0.2 2",$1 }' > param_R.dat
ls *EHZ.sac *EHR.sac | awk '{print "0 2.5 5 10 100 20 1 0.5 0.2 2",$1 }' > param_R.dat

./aftani_c_pgl_output_amp.Utah/aftani_c_pgl_test param_R.dat 
                                                                                           
#lf_spec_snr_rms_fast_v4 sac_R.lst 0.1 30 1 4.5                                                                                                                          
#lf_spec_snr_rms_fast_v4 sac_R.lst 10 100 2.5 5                                                                 


#cd ..
#end

