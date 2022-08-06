%% New test% reconstruct phi from non-rotation free mesh
clc
clear
close all
format long
%% Uniform grid
%%
N=97;m=N;n=N;h1=2; hh=2; showM = 512;
y=1:m;
[X,Y]=ndgrid(y,y);
X_new=X;Y_new=Y;
Blank_white=zeros(N)+255;
%% construct transformations

load Phi1_circ
load Phi2_circ
load cv_circ
load jd_circ

phi1=Phi1_circ;
phi2=Phi2_circ;
jd_phi=jd_circ;
cv_phi=cv_circ;

load Phi1_tbl
load Phi2_tbl
load cv_tbl
load jd_tbl

psi1=Phi1_tbl;
psi2=Phi2_tbl;
jd_psi=jd_tbl;
cv_psi=cv_tbl;

load Phi1_tile
load Phi2_tile
load cv_tile
load jd_tile

psi1_trans=Phi1_tile;
psi2_trans=Phi2_tile;
jd_psi_trans=jd_tile;
cv_psi_trans=cv_tile;

h=2;
figure(1)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-phi1(1:h:end,i),phi2(1:h:end,i),'k-'), hold on
    plot(-phi1(i,1:h:end),phi2(i,1:h:end),'k-'), hold on
end
axis([-N,-1,1,N]);
figure(2)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-psi1(1:h:end,i),psi2(1:h:end,i),'k-'), hold on
    plot(-psi1(i,1:h:end),psi2(i,1:h:end),'k-'), hold on
end
axis([-N,-1,1,N]); 
figure(3)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-psi1_trans(1:h:end,i),psi2_trans(1:h:end,i),'k-'), hold on
    plot(-psi1_trans(i,1:h:end),psi2_trans(i,1:h:end),'k-'), hold on
end
axis([-N,-1,1,N]); 
%% fwd
tic;
[phi1_fw,phi2_fw,U1_fw,U2_fw,phi_cur1_fw, phi_cur2_fw,U1_cur_fw,U2_cur_fw]=PJDC_on_given_mesh2(jd_psi_trans,cv_psi_trans,N,phi1,phi2);
toc;

[JD_phi_fw, CV_phi_fw]=compute_JD_and_Curl(phi1_fw,phi2_fw,1);
diff_jd_phi_fw=abs(JD_phi_fw-jd_psi_trans);
diff_cv_phi_fw=abs(CV_phi_fw-cv_psi_trans);
diff_jd_phi_max_fw=max(max(diff_jd_phi_fw))
diff_cv_phi_max_fw=max(max(diff_cv_phi_fw))

diff_Tx_phi_fw=phi1_fw-psi1_trans;
diff_Ty_phi_fw=phi2_fw-psi2_trans;
diff_mag_phi_fw=(diff_Tx_phi_fw.^2+diff_Ty_phi_fw.^2).^(0.5);
diff_max_mag_phi_fw=max(max(diff_mag_phi_fw))
rel_diff_max_mag1_phi_fw=max(max(diff_mag_phi_fw./((psi1_trans.^2+psi1_trans.^2).^(0.5))))
rel_diff_max_mag2_phi_fw=max(max(diff_mag_phi_fw./(100^2)))

figure(4)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-phi_cur1_fw(1:h:end,i),phi_cur2_fw(1:h:end,i),'r-'), hold on
    plot(-phi_cur1_fw(i,1:h:end),phi_cur2_fw(i,1:h:end),'r-'), hold on
end
axis([-N,-1,1,N]);
figure(5)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
quiver(-X(1:1:N, 1:1:N), Y(1:1:N, 1:1:N), -U1_cur_fw(1:1:N, 1:1:N), U2_cur_fw(1:1:N, 1:1:N), 0);
axis([-N,-1,1,N]);  
figure(6)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-psi1_trans(1:h:end,i), psi2_trans(1:h:end,i),'k-'), hold on
    plot(-psi1_trans(i,1:h:end), psi2_trans(i,1:h:end),'k-'), hold on
end
for i = 1:h:N
    plot(-phi1_fw(1:h:end,i), phi2_fw(1:h:end,i),'r-'), hold on
    plot(-phi1_fw(i,1:h:end), phi2_fw(i,1:h:end),'r-'), hold on
end
axis([-N,-1,1,N]); 
%% transitivity %% fwd 
tic;
[phi1_fw2,phi2_fw2,U1_fw2,U2_fw2,phi_cur1_fw2, phi_cur2_fw2,U1_cur_fw2,U2_cur_fw2]=PJDC_on_given_mesh2(jd_psi,cv_psi,N,psi1_trans,psi2_trans);
toc;

[JD_phi_fw2, CV_phi_fw2]=compute_JD_and_Curl(phi1_fw2,phi2_fw2,1);
diff_jd_phi_fw2=abs(JD_phi_fw2-jd_psi);
diff_cv_phi_fw2=abs(CV_phi_fw2-cv_psi);
diff_jd_phi_max_fw2=max(max(diff_jd_phi_fw2))
diff_cv_phi_max_fw2=max(max(diff_cv_phi_fw2))

diff_Tx_phi_fw2=phi1_fw2-psi1;
diff_Ty_phi_fw2=phi2_fw2-psi2;
diff_mag_phi_fw2=(diff_Tx_phi_fw2.^2+diff_Ty_phi_fw2.^2).^(0.5);
diff_max_mag_phi_fw2=max(max(diff_mag_phi_fw2))
rel_diff_max_mag1_phi_fw2=max(max(diff_mag_phi_fw2./((psi1.^2+psi2.^2).^(0.5))))
rel_diff_max_mag2_phi_fw2=max(max(diff_mag_phi_fw2./(100^2)))

figure(7)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot(-phi_cur1_fw2(i,1:hh:N), phi_cur2_fw2(i,1:hh:N),'r-'),hold on
    plot(-phi_cur1_fw2(1:hh:N,i), phi_cur2_fw2(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 

figure(8)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
quiver(-X(1:1:N, 1:1:N), Y(1:1:N, 1:1:N), -U1_cur_fw2(1:1:N, 1:1:N), U2_cur_fw2(1:1:N, 1:1:N), 0);
axis([-N,-1,1,N]); 
figure(9)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-psi1(1:h:end,i),psi2(1:h:end,i),'k-'), hold on
    plot(-psi1(i,1:h:end),psi2(i,1:h:end),'k-'), hold on
end
for i = 1:hh:N
    plot(-phi1_fw2(i,1:hh:N), phi2_fw2(i,1:hh:N),'r-'),hold on
    plot(-phi1_fw2(1:hh:N,i), phi2_fw2(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 


transt_check_x = interpn(X, Y, phi_cur1_fw, phi_cur1_fw2, phi_cur2_fw2,'makima'); 
transt_check_y = interpn(X, Y, phi_cur2_fw, phi_cur1_fw2, phi_cur2_fw2,'makima'); 

transt_check_phi1 = interpn(X, Y, phi1, transt_check_x, transt_check_y,'makima'); 
transt_check_phi2 = interpn(X, Y, phi2, transt_check_x, transt_check_y,'makima'); 


figure(10)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot(-transt_check_x(i,1:hh:end), transt_check_y(i,1:hh:end),'r-'),hold on
    plot(-transt_check_x(1:hh:end,i), transt_check_y(1:hh:end,i),'r-'),hold on
end
axis([-N,-1,1,N]); 
figure(11)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot(-transt_check_phi1(i,1:hh:N), transt_check_phi2(i,1:hh:N),'r-'),hold on
    plot(-transt_check_phi1(1:hh:N,i), transt_check_phi2(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 

%% jump forward
tic;
[phi1_jfw,phi2_jfw,U1_jfw,U2_jfw,phi_cur1_jfw, phi_cur2_jfw,U1_cur_jfw,U2_cur_jfw]=PJDC_on_given_mesh2(jd_psi,cv_psi,N,phi1,phi2);
toc;

[JD_phi_jfw, CV_phi_jfw]=compute_JD_and_Curl(phi1_jfw,phi2_jfw,1);
diff_jd_phi_jfw=abs(JD_phi_jfw-jd_psi);
diff_cv_phi_jfw=abs(CV_phi_jfw-cv_psi);
diff_jd_phi_max_jfw=max(max(diff_jd_phi_jfw))
diff_cv_phi_max_jfw=max(max(diff_cv_phi_jfw))

diff_Tx_phi_jfw=phi1_jfw-psi1;
diff_Ty_phi_jfw=phi2_jfw-psi2;
diff_mag_phi_jfw=(diff_Tx_phi_jfw.^2+diff_Ty_phi_jfw.^2).^(0.5);
diff_max_mag_phi_jfw=max(max(diff_mag_phi_jfw))
rel_diff_max_mag1_phi_jfw=max(max(diff_mag_phi_jfw./((psi1.^2+psi2.^2).^(0.5))))
rel_diff_max_mag2_phi_jfw=max(max(diff_mag_phi_jfw./(100^2)))

figure(12)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot( -phi_cur1_jfw(i,1:hh:N), phi_cur2_jfw(i,1:hh:N),'r-'),hold on
    plot( -phi_cur1_jfw(1:hh:N,i), phi_cur2_jfw(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 

figure(13)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
quiver(-X(1:1:N, 1:1:N), Y(1:1:N, 1:1:N), -U1_cur_jfw(1:1:N, 1:1:N), U2_cur_jfw(1:1:N, 1:1:N), 0);
axis([-N,-1,1,N]); 

figure(14)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:h:N
    plot(-psi1(1:h:end,i),psi2(1:h:end,i),'k-'), hold on
    plot(-psi1(i,1:h:end),psi2(i,1:h:end),'k-'), hold on
end
for i = 1:hh:N
    plot(-phi1_jfw(i,1:hh:N), phi2_jfw(i,1:hh:N),'r-'),hold on
    plot(-phi1_jfw(1:hh:N,i), phi2_jfw(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 

figure(15)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot(-transt_check_x(i,1:hh:end), transt_check_y(i,1:hh:end),'k-'),hold on
    plot(-transt_check_x(1:hh:end,i), transt_check_y(1:hh:end,i),'k-'),hold on
end
for i = 1:hh:N
    plot( -phi_cur1_jfw(i,1:hh:N), phi_cur2_jfw(i,1:hh:N),'r-'),hold on
    plot( -phi_cur1_jfw(1:hh:N,i), phi_cur2_jfw(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 
figure(16)
imshow(imresize(Blank_white, [showM showM]), [], 'border', 'tight'), hold on
for i = 1:hh:N
    plot(-transt_check_phi1(i,1:hh:N), transt_check_phi2(i,1:hh:N),'k-'),hold on
    plot(-transt_check_phi1(1:hh:N,i), transt_check_phi2(1:hh:N,i),'k-'),hold on
end
for i = 1:hh:N
    plot(-phi1_jfw(i,1:hh:N), phi2_jfw(i,1:hh:N),'r-'),hold on
    plot(-phi1_jfw(1:hh:N,i), phi2_jfw(1:hh:N,i),'r-'),hold on
end
axis([-N,-1,1,N]); 
%  ts: 9.9938e-11 r: 0.00128 ei: 17367 ti: 43921
% Elapsed time is 421.593143 seconds.

% diff_jd_phi_max_fw =
% 
%    0.033277909902358
% 
% 
% diff_cv_phi_max_fw =
% 
%    0.022229759769765
% 
% 
% diff_max_mag_phi_fw =
% 
%    0.250763947308684
% 
% 
% rel_diff_max_mag1_phi_fw =
% 
%    0.069241365252727
% 
% 
% rel_diff_max_mag2_phi_fw =
% 
%      2.507639473086840e-05
% 
%  ts: 9.9962e-11 r: 0.0019999 ei: 21126 ti: 51435
% Elapsed time is 470.314256 seconds.
% 
% diff_jd_phi_max_fw2 =
% 
%    0.038555819391753
% 
% 
% diff_cv_phi_max_fw2 =
% 
%    0.032688657031473
% 
% 
% diff_max_mag_phi_fw2 =
% 
%    0.341425569049768
% 
% 
% rel_diff_max_mag1_phi_fw2 =
% 
%    0.007680656137979
% 
% 
% rel_diff_max_mag2_phi_fw2 =
% 
%      3.414255690497683e-05
% 
%  ts: 9.9919e-11 r: 0.00096487 ei: 20559 ti: 50302
% Elapsed time is 462.600204 seconds.
% 
% diff_jd_phi_max_jfw =
% 
%    0.028323414200295
% 
% 
% diff_cv_phi_max_jfw =
% 
%    0.017404427429579
% 
% 
% diff_max_mag_phi_jfw =
% 
%    0.144427922944189
% 
% 
% rel_diff_max_mag1_phi_jfw =
% 
%    0.005195009468293
% 
% 
% rel_diff_max_mag2_phi_jfw =
% 
%      1.444279229441889e-05
% 

