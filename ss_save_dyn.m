function ss_save_dyn(g,rf,ang,thk, isodelay, format, fspec, a_angs, root_fname, pulse_num)
% SS_SAVE - Save spectral-spatial pulse 
% Uses Chuck Cunningham's format for GE systems, and creates associated
% .dat-file
% Pulse parameters saved in header for Varian and Bruker files
%
%  ss_save_dyn(g,rf,ang,thk, format, fspec, a_angs, root_fname, pulse_num)
%
%  g - in G/cm    
%  rf - in G
%  ang - flip angle in radians
%  thk - thickness in cm
%  isodelay (optional) - delay from in-phase point to end of pulse (GE
%  definition)
%  format (optional) - 'GE' (default), 'Varian', 'Bruker'
%  fspec (optional) - frequency bands (Hz) to write in file
%  a_angs (optional) - band amplidutes (radians) to write in file
%  root_fname - root file name (no prompting)
%  pulse_num - pulse number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2011 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss_globals;
    
if (nargin < 5) || isempty(isodelay)
    isodelay = length(rf)*SS_TS/ 2;
end

if (nargin < 6) || isempty(format)
    format = 'GE';
end

   switch format,
        case 'GE'
            % force even number of samples due to some GE sequences having
            % issues loading odd number of samples
            if rem(length(rf),2)
                rf = [rf(:); 0];
                g = [g(:); 0];
            end
            
        case 'Varian'
        case 'Bruker'
        otherwise
            error(sprintf(['Format save type of: %s not' ...
                ' recognized'], format));
    end;
    
    maxg = max(abs(g));
    if maxg ~= 0, 
	gn = g/maxg;
    else
	gn = g;
    end;
    
    maxrf = max(abs(rf));
    rfn = rf / maxrf;

    
    % calculate pulse parameters
    nrf = length(rf);
    abswidth = sum(abs(rfn))/nrf;
    effwidth = sum(abs(rfn).^2)/nrf;
    area = sum(abs(rfn))/nrf;
    pon = (rfn >= 0.00001);
    temp_pw = 0;
    max_pw = 0;
    for n=1:nrf
	temp_pw = temp_pw + pon(n);
	if (and(pon(n) == 0, temp_pw ~= 0))
	    max_pw = max(max_pw, temp_pw);
	    temp_pw = 0;
	end;
    end;
    max_pw = max_pw / n;
    
    dty_cyc = sum(abs(rfn) > 0.2236)/nrf; 
    if dty_cyc < max_pw, 
	dty_cyc = max_pw;
    end;

    max_b1 = max(abs(rf));
    int_b1_sqr = sum(abs(rf).^2 * SS_TS * 1e3);
    rms_b1 = sqrt(sum(abs(rf).^2))/nrf;
    thk_scale = thk * SS_GAMMA / SS_GAMMA_HYDROGEN * 10;
    
    
    % Allow magnitude of RF to go negative, this will
    % help reduce sensitivity to theta modulation since
    % there is a possible delay on the system for this case
    %
    
    % Find pi jumps in phase and remove
    % Do this by doubling phase and unwrapping 2*pi jumps
    %
    dang_rf = 2*angle(rfn);
    dang_rf = dang_rf - dang_rf(1);
    dang_rf = unwrap(dang_rf, 0.98*pi);
    ang_rf = mod(dang_rf/2+pi,2*pi)-pi;
    
    mag_rf = real(rfn .* exp(-i*ang_rf));
    
%     if 1, 
% 	figure;
% 	subplot(211)
% 	t = 1:length(rf);
% 	plot(abs(rf));
% 	hold on;
% 	plot(maxrf*mag_rf,'r--')
% 	
% 	subplot(212);
% 	plot(real(rf));
% 	hold on;
% 	plot(imag(rf), 'b--');
% 	plot(real(maxrf*mag_rf.*exp(i*ang_rf)), 'r:');
% 	plot(imag(maxrf*mag_rf.*exp(i*ang_rf)), 'r:');
%     end;

    
    
    switch (format)
        case 'GE'
    
            dat_name = sprintf('%s%03d.dat', root_fname, pulse_num);
            fid = fopen(dat_name, 'w');
            if fid == -1, 
            fprintf(1, 'Error opening %s \n', dat_name);
            return;
            end;

            fprintf(fid,'%10d \t\t #Spectral-Spatial\n', 1);
            fprintf(fid,'%10d \t\t #res\n', length(rf));
            fprintf(fid,'%10d \t\t #pw\n',round(length(rf)*SS_TS*1e6));
            fprintf(fid,'%10.7f \t\t #nom_flip \n',ang*180/pi);
            fprintf(fid,'%10.7f \t\t #abswidth \n',abswidth);
            fprintf(fid,'%10.7f \t\t #effwidth \n',effwidth);
            fprintf(fid,'%10.7f \t\t #area \n',area);
            fprintf(fid,'%10.7f \t\t #dtycyc \n',dty_cyc);
            fprintf(fid,'%10.7f \t\t #maxpw \n',max_pw);
            gamscale = SS_GAMMA/SS_GAMMA_HYDROGEN; % GE assumes max B1 is for application at 1H
            fprintf(fid,'%10.7f \t\t #max_b1 \n',max_b1 * gamscale); 
            fprintf(fid,'%10.7f \t\t #max_int_b1_sqr \n',int_b1_sqr* gamscale^2);
            fprintf(fid,'%10.7f \t\t #max_rms_b1 \n',rms_b1* gamscale^2);
            fprintf(fid,'%10.3f \t\t #a_gzs \n',maxg);
            fprintf(fid,'%10.3f \t\t #nom_thk(mm) \n',thk * gamscale * 10);
            fprintf(fid,'%10d \t\t #isodelay\n',round(isodelay*1e6));
            fprintf(fid,'%10d \t\t #g_pow \n',0);
            fprintf(fid,'%10d \t\t #g_pos_pow \n',0);
            fprintf(fid,'%10d \t\t #g_neg_pow \n',0);
            fprintf(fid,'%10d \t\t #g_abs \n',0);
            fprintf(fid,'%10d \t\t #g_dgdt \n',0);
            fprintf(fid,'%10d \t\t #g_pwm \n',0);
            fprintf(fid,'%10d \t\t #g_pwm_abs \n',0);
            fprintf(fid,'# *************************************\n');
            if (nargin  > 6)
                for b = 1:length(a_angs)
                   fprintf(fid,'# Band %d: [%.2f, %.2f] Hz, %.2f degree flip\n', ...
                       b, fspec(2*b-1), fspec(2*b), a_angs(b)*180/pi);
                end
            end
            
            fclose(fid);

            % Now write out RF and Gradient 
            %
            rho_fname = sprintf('%s%03d.rho', root_fname, pulse_num);
            signa(mag_rf,rho_fname,1);

            theta_fname = sprintf('%s%03d.pha', root_fname, pulse_num);
            signa(-ang_rf,theta_fname,1/pi);

            g_fname = sprintf('%s%03d.grd', root_fname, pulse_num);
            signa(gn,g_fname,1);
    
        case 'Varian'
            
            fid = fopen(sprintf('%s%03d.RF', root_fname, pulse_num),'wt');
            fprintf(fid,'# %s\n', sprintf('%s%03d.RF', root_fname, pulse_num));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# Spectral-spatial Matlab Package\n');
            fprintf(fid,'# Nucleus = %s\n',SS_NUCLEUS);
            fprintf(fid,'# Duration = %d us\n',round(length(rf)*SS_TS*1e6));
            fprintf(fid,'# Isodelay = %d us\n',round(isodelay*1e6));
            fprintf(fid,'# Resolution = %d us\n', SS_TS*1e6);
            fprintf(fid,'# Flip = %.2f degrees\n',ang*180/pi);
            fprintf(fid,'# Max B1 = %.4f Gauss\n',max_b1);
             if (nargin > 6)
                for b = 1:length(a_angs)
                   fprintf(fid,'# Band %d: [%.2f, %.2f] Hz, %.2f degree flip\n', ...
                       b, fspec(2*b-1), fspec(2*b), a_angs(b)*180/pi);
                end
            end
           fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# VERSION       Matlab\n');
            fprintf(fid,'# TYPE          selective\n');
            fprintf(fid,'# MODULATION    amplitude\n');
            fprintf(fid,'# EXCITEWIDTH   -1.0000\n');
            fprintf(fid,'# INVERTWIDTH   -1.0000\n');  
            fprintf(fid,'# INTEGRAL      %1.4f\n',    sum(abs(rfn))/length(rfn)); 
            fprintf(fid,'# RF_FRACTION   -1.0000\n');
            fprintf(fid,'# STEPS         %d\n',length(rfn));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'%3.2f  %4.2f  1.0\n',[-angle(rfn(:)).'*180/pi; 1023*abs(rfn(:)).']);    
            fclose(fid);

            fid = fopen(sprintf('%s%03d.GRD', root_fname, pulse_num),'wt');
            fprintf(fid,'# %s\n', sprintf('%s.GRD', root_fname));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# Spectral-spatial Matlab Package\n');
            fprintf(fid,'# Nucleus = %s\n',SS_NUCLEUS);
            fprintf(fid,'# Duration = %d us\n',round(length(g)*SS_TS*1e6));
            fprintf(fid,'# Resolution = %d us\n', SS_TS*1e6);
            fprintf(fid,'# Points = %d\n',length(g));
            fprintf(fid,'# Max Gradient Strength = %.4f Gauss/cm\n',maxg);
            fprintf(fid,'# (Max Gradient Strength Constraint = %.2f Gauss/cm)\n',SS_MXG);
            fprintf(fid,'# Max Slew Rate = %.2f Gauss/cm/ms\n',SS_MXS);
            fprintf(fid,'# Slice thickness = %.1f mm\n',thk*10);
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'%d  1\n',round(32767*gn));
            fclose(fid);
    


        case 'Bruker'
            
            fid = fopen(sprintf('%s%03d.exc', root_fname, pulse_num),'wt');
  
            fprintf(fid,'##TITLE= %s\n', sprintf('%s%03d.RF', root_fname, pulse_num));
            fprintf(fid,'##JCAMP-DX= 5.00 BRUKER JCAMP library\n');
            fprintf(fid,'##DATA TYPE= Shape Data\n');
            fprintf(fid,'##ORIGIN= BRUKER MEDICAL\n');
            fprintf(fid,'##OWNER= <BRUKER MEDICAL>\n');
            fprintf(fid,'##DATE= 20171106\n');
            fprintf(fid,'##MINX= \n');
            fprintf(fid,'##MAXX= \n');
            fprintf(fid,'##MINY= \n');
            fprintf(fid,'##MAXY= \n');
            fprintf(fid,'##$SHAPE_EXMODE=  Excitation\n');
            fprintf(fid,'##$SHAPE_TOTROT= 90.000000e+00\n');
            fprintf(fid,'##$SHAPE_BWFAC= %6.6f\n',z_tb);
            fprintf(fid,'##$SHAPE_INTEGFAC=      %1.4f\n',  sum(cos(-angle(rfn)).*abs(rfn))/length(rfn)); %changed by Vickie for Varian 

            fprintf(fid,'##$SHAPE_REPHFAC= 50\n');
            fprintf(fid,'##$SHAPE_TYPE=conventional\n');
            fprintf(fid,'##$SHAPE_MODE= 0 \n');
            fprintf(fid,'##NPOINTS= %d\n',length(g)); 
            
            fprintf(fid,'##Spectral_spatial_Matlab_Package\n');
            fprintf(fid,'##Nucleus= %s\n',SS_NUCLEUS);
            fprintf(fid,'##Duration= %d us\n',round(length(rf)*SS_TS*1e6));
            fprintf(fid,'##Resolution= %d us\n', SS_TS*1e6);
            fprintf(fid,'##Flip= %.2f degrees\n',ang*180/pi);
            fprintf(fid,'##Max_B1 = %.4f Gauss\n',max_b1);
            fprintf(fid,'##Spatial_time_bandwidth= %d\n',z_tb);
            if SS_SLR_FLAG == 0
               fprintf(fid,'##SHAPE= small tip\n');
            else
               fprintf(fid,'##SHAPE= SLR\n');
            end
             if (nargin > 6)
                for b = 1:length(a_angs)
                   fprintf(fid,'##Band_%d= [%.2f, %.2f] Hz, %.2f degree flip, %.3f ripple\n', ...
                       b, fspec(2*b-1), fspec(2*b), a_angs(b)*180/pi, d(b)/sin(max(a_angs)));
                end
            end             
            fprintf(fid,'##XYPOINTS= (XY..XY)\n');
            
            %fprintf(fid,'%3.2f  %4.2f  1.0\n',[-angle(rfn(:)).'*180/pi; 1023*abs(rfn(:)).']);    
            %fprintf(fid,'%5.2f  %5.2f \n',[-angle(rfn(:)).'*180/pi; 1023*abs(rfn(:)).']);    
            fprintf(fid,'%5.2f, %5.2f \n',[abs(rfn(:)).'*100'; -angle(rfn(:)).'*180/pi' ]);    
       

            fclose(fid);

            fid = fopen(sprintf('%s%03d.gp', root_fname, pulse_num),'wt');
            fprintf(fid,'##TITLE= %s\n', sprintf('%s%03d.gp', root_fname, pulse_num));
            fprintf(fid,'## Spectral_spatial_Matlab_Package\n');
            fprintf(fid,'##Nucleus= %s\n',SS_NUCLEUS);
            fprintf(fid,'##SHAPE= %s\n',ss_type);
            fprintf(fid,'##Duration= %d us\n',round(length(g)*SS_TS*1e6));
            fprintf(fid,'##Resolution= %d us\n', SS_TS*1e6);
            fprintf(fid,'##Points= %d\n',length(g));
            fprintf(fid,'##Max_Gradient_Strength= %.4f Gauss/cm\n',maxg);
            fprintf(fid,'##(Max_Gradient_Strength_Constraint= %.2f Gauss/cm)\n',SS_MXG);
            fprintf(fid,'##Max_Slew Rate= %.2f Gauss/cm/ms\n',SS_MXS);
            fprintf(fid,'##Slice_thickness= %.1f mm\n',thk*10);
            fprintf(fid,'%5.3f  \n',g);
            fclose(fid);
    end
    
    % Data Acquisition Descriptor (DAD) XML file containing RF pulse information
    % for more info see SIVIC project https://github.com/SIVICLab/sivic
    dad_name = sprintf('%s.xml', root_fname);
    
    if pulse_num ==1
        docNode = com.mathworks.xml.XMLUtils.createDocument('svk_data_acquisition_description');
        dad = docNode.getDocumentElement;
        
        dad_version = docNode.createElement('version');
        dad_version.appendChild(docNode.createTextNode('0'));
        dad.appendChild(dad_version);
        
        encoding =docNode.createElement('encoding');
        dad.appendChild(encoding);
        
        excitation =docNode.createElement('excitation');
        encoding.appendChild(excitation);
        
        spectralType =docNode.createElement('spectralType');
        spectralType.appendChild(docNode.createTextNode('selective'));
        excitation.appendChild(spectralType);
        
        spatialType =docNode.createElement('spatialType');
        spatialType.appendChild(docNode.createTextNode('selective'));
        excitation.appendChild(spatialType);
        
        pulseName =docNode.createElement('pulseName');
        pulseName.appendChild(docNode.createTextNode(root_fname));
        excitation.appendChild(pulseName);
        
        % flip angle(s) and associated frequency bands
        if (nargin > 6)
            for b = 1:length(a_angs)
                curr_node = docNode.createElement('flipAngle_deg');
                
                % need to convert to strings??
                curr_node.setAttribute('frequencyMin_Hz',num2str(min(fspec(2*b-1),fspec(2*b))));
                curr_node.setAttribute('frequencyMax_Hz',num2str(max(fspec(2*b-1),fspec(2*b))));
                
                curr_node.appendChild(docNode.createTextNode(num2str(a_angs(b)*180/pi)));
                
                excitation.appendChild(curr_node);
            end
        end
        
    else
        docNode = xmlread(dad_name);
        
        % flip angle(s) and associated frequency bands
        if (nargin > 6)
            flipAngle_deg = docNode.getElementsByTagName('flipAngle_deg');
            
            for b = 1:length(a_angs)
                curr_node = flipAngle_deg.item(b-1);
                curr_node_child = curr_node.getFirstChild;
                new_flipAngle_deg_data = [char(curr_node_child.getData), ' ' num2str(a_angs(b)*180/pi)];
                curr_node.replaceChild(docNode.createTextNode(new_flipAngle_deg_data), curr_node_child);
            end
        end
        
    end
    
    
    % add pulse frequency
    
    xmlwrite(dad_name,docNode);
    
    
    
    
    

    
