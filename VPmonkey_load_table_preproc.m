

% load table with preprocessing info
function tbl_preproc = VPmonkey_load_table_preproc

    temp = readtable([getRoot filesep 'VPmonkey/preprocessing_by_session.xls'], ...
    "UseExcel", false);
    
    tbl_preproc = temp(1:end,[1,2,14:16 22]);
    tbl_preproc.Properties.VariableNames = {'sesh','sub','rej_aud','rej_som','rej_vis','rej_eye'};

end