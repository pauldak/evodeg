% Generates a PyMOL batch script that globally aligns all proteins to the reference organism
% (S. cerevisiae) and saves transformed coordinates (_tx.pdb) and alignment files (_tx.aln).

function pdbtx_all_protein_initation()
    load('folder_paths.mat');
    number_iterations = size(data_table, 1);
    %creating pymol file
    pymol_script_filename = 'create_pdbtx.pml';
    fid = fopen(pymol_script_filename, 'w');
    for i = 1: number_iterations
        crnt_folder = [data_folder data_table.Gene_Group{i} '\'];
        crnt_zero_pdb_path = [crnt_folder zero_organism '\' data_table.(zero_organism){i} '.pdb'];
        for k = 1: length(organism_list)           
            crnt_organism = organism_list{k};
            crnt_organism_folder = [crnt_folder crnt_organism '\'];
            crnt_prot_id = data_table.(crnt_organism){i};
            crnt_org_pdb_path = [crnt_organism_folder crnt_prot_id '.pdb'];
            %writing comands to pml file
            fprintf(fid, 'load %s, zero_prot\n', crnt_zero_pdb_path);
            fprintf(fid, 'load %s, crnt_prot\n', crnt_org_pdb_path);
            fprintf(fid, 'align crnt_prot, zero_prot, object="crnt_align"\n');
            fprintf(fid, 'zoom complete=1\n');
            fprintf(fid, 'save %s%s_tx.pdb, crnt_prot\n', crnt_organism_folder, crnt_prot_id);
            fprintf(fid, 'save %s%s_tx.aln, crnt_align\n', crnt_organism_folder, crnt_prot_id);
            fprintf(fid, 'png %s%s_tx, dpi=900\n', crnt_organism_folder, crnt_prot_id);
            fprintf(fid, 'save %s%s_tx.pse\n', crnt_organism_folder, crnt_prot_id);
            fprintf(fid, 'reinitialize\n');
        end
    end
    fclose(fid);
    pause(5); %let the file be saved
    pymol_adrs = '"C:\Users\SEA\AppData\Local\Schrodinger\PyMOL2\PyMOLWin.exe"';
    % Execute the PyMOL script
    pymol_command = [pymol_adrs ' -cq ' pymol_script_filename];
    status = system(pymol_command);
    if status == 0
        disp('PyMOL opened with the PML script.');
    else
        error('Failed to open PyMOL with the PML script.');
    end
end


