function [Mw, pISO, pDC, pCLVD, Bplunge]=RadPatt(input, input_flag, plot_flag)
  % Simple program to decompose and display the radiation pattern for a moment tensor.
  %
  % For conventions see Aki & Richards (2002).
  % e.g., X - North, Y - East, Z - Down.
  %
  % For P-waves red indicates outward radiation, compressional first P motion, and is black on a beachball.
  % For P-waves blue indicates inward radiation, dilitational first P motion, and is white on a beachball.
  % For Sh-waves we use right-handed coordinates. Red indicates clock-wise shearing, when looking into the ground and vice versa.
  % P-axis depicted by white diamonds.
  % B-axis depicted by black x's.
  % T-axis depicted by black circles.
  %
  %
  % Strike angle is measured clock-wise from north and defined as the fault intersection with surface (0-360 degrees).
  % Strike is defined such that the fault dips to the right when facing the strike direction, i.e., the hanging wall is on right side.
  % Dip angle is measure as angle from horizontal down to the fault plane (0-90 degrees).
  % Rake(slip) angle is the direction the hanging wall moves on the fault plane during slip measured from horizontal (-180 to 180 degrees).
  % A rake of 0 degrees means the hanging wall moved in the strike direction (left-lateral/sinistral motion), and vice versa.
  % When rake is greater than zero the hanging wall moved up (thrust/reverse motion), and vice versa.
  %
  % References:
  % 
  % Aki, K., & Richards, P. G. (2002). Quantitative Seismology.
  % Frohlich, C. (1992). Triangle diagrams: ternary graphs to display similarity and diversity of earthquake focal mechanisms. Physics of the Earth and Planetary Interiors, 75(1), 193-198, doi: 10.1016/0031-9201(92)90130-N
  % Frohlich, C. (1996). Cliff's nodes concerning plotting nodal lines for P, Sh and Sv. Seismological Research Letters, 67(1), 16-24, doi: 10.1785/gssrl.67.1.16.
  % Hudson, J. A., Pearce, R. G., & Rogers, R. M. (1989). Source type plot for inversion of the moment tensor. Journal of Geophysical Research: Solid Earth (1978?2012), 94(B1), 765-774, doi: 10.1029/JB094iB01p00765.
  % Jost, M. U., & Herrmann, R. B. (1989). A student's guide to and review of moment tensors. Seismological Research Letters, 60(2), 37-57, doi: 10.1785/gssrl.60.2.37.
  % Vavrycuk, V. (2015). Moment tensor decomposition revisited. Journal of Seismology, 19.1, 231-252, doi: 10.1007/s10950-014-9463-y.
  %
  % Written by Ryan Schultz.
  
  %%% STILL HAS A BUG where orientations along Z axis don't work. or maybe
  %%% that's actually the way things work...?
  
  
  % Initialize variables destined for output.
  Mw=-99;
  pISO=-99;
  pDC=-99;
  pCLVD=-99;
  
  % Check input makes sense
  if(max(isnan(input)))
      fprintf('Improper input fields (i.e., Inf or NaN).  Aborting.\n');
      return;
  elseif(max(isinf(input)))
      fprintf('Improper input fields (i.e., Inf or NaN).  Aborting.\n');
      return;
  end;
  
  % Construct moment tensor based on flagged input.
  if(input_flag==1) % Cartesian inputs.
      
      % Check that input length is proper.
      if(length(input)~=6)
          fprintf('Improper input length for Cartesian tensor elements.  Aborting.\n');
          return;
      end;
      Mxx=input(1);
      Mxy=input(2);
      Mxz=input(3);
      Myy=input(4);
      Myz=input(5);
      Mzz=input(6);
      
      % GMT CMT formatting
      Mrr=Mzz;
      Mrdelta=Mxz;
      Mrphi=-Myz;
      Mdeltadelta=Mxx;
      Mdeltaphi=-Mxy;
      Mphiphi=Myy;
      fprintf('%e %e %e %e %e %e\n', Mrr, Mdeltadelta, Mphiphi, Mrdelta, Mrphi, Mdeltaphi);
      
  elseif(input_flag==2) % Strike, dip, rake inputs.
      
      % Check that input length is proper.
      if(length(input)~=3)
          fprintf('Improper input length for strike,dip, and rake.  Aborting.\n');
          return;
      end;
      strike=input(1);
      dip=input(2);
      rake=input(3);
      
      % Convert to Cartesian tensor elements (Aki & Richards (2002), p 112-113).
      Mxx=-( sind(dip)*cosd(rake)*sind(2*strike)  +     sind(2*dip)*sind(rake)*(sind(strike)^2) );
      Mxy=+( sind(dip)*cosd(rake)*cosd(2*strike)  + 0.5*sind(2*dip)*sind(rake)*sind(2*strike)   );
      %Myx=Mxy;
      Mxz=-( cosd(dip)*cosd(rake)*cosd(strike)    +     cosd(2*dip)*sind(rake)*sind(strike)     );
      %Mzx=Mxz;
      Myy=+( sind(dip)*cosd(rake)*sind(2*strike)  -     sind(2*dip)*sind(rake)*(cosd(strike)^2) );
      Myz=-( cosd(dip)*cosd(rake)*sind(strike)    -     cosd(2*dip)*sind(rake)*cosd(strike)     );
      %Mzy=Myz;
      Mzz=sind(2*dip)*sind(rake);
      
  elseif(input_flag==3) % Spherical inputs.
      
      % Check that input length is proper.
      if(length(input)~=6)
          fprintf('Improper input length for Spherical tensor elements.  Aborting.\n');
          return;
      end;
      Mrr=input(1);
      Mrdelta=input(2);
      Mrphi=input(3);
      Mdeltadelta=input(4);
      Mdeltaphi=input(5);
      Mphiphi=input(6);
      
      % Convert to spherical tensor elements (Aki & Richards (2002), p 112-113).
      %Mrr=Mzz;
      %Mrdelta=Mxz;
      %Mrphi=-Myz;
      %Mdeltadelta=Mxx;
      %Mdeltaphi=-Mxy;
      %Mphiphi=Myy;
      
      % Convert to Cartesian tensor elements (Aki & Richards (2002), p 112-113).
      Mzz=Mrr;
      Mxz=Mrdelta;
      Myz=-Mrphi;
      Mxx=Mdeltadelta;
      Mxy=-Mdeltaphi;
      Myy=Mphiphi;
      
      % GMT CMT formatting
      fprintf('%e %e %e %e %e %e\n', Mrr, Mdeltadelta, Mphiphi, Mrdelta, Mrphi, Mdeltaphi);
      
  else
      fprintf('Improper input flag.  Aborting.\n');
      return;
  end
  
  % Put tensor elements into a tensor.
  M=[Mxx, Mxy, Mxz; Mxy, Myy, Myz; Mxz, Myz, Mzz];
  
  % Compute seismic moment.
  Mo=sqrt(sum(sum(M.^2))/2);
  Mw=(2/3)*(log10(Mo)-16.1);
  
  % Diagonalization of moment tensor.
  [V,D]=eig(M);
  lambda=diag(D);
  
  % Re-sort the order of eigenvalues and eigenvectors.
  [m,I]=sort(lambda,'descend');
  V=[V(:,I(1)),V(:,I(2)),V(:,I(3))];
  D=diag(m);
  lambda=[lambda(I(1)); lambda(I(2)); lambda(I(3))];
  c=max(abs(lambda));  clim=[-c c];
  
  % Compute relative proportions of isotropic, double couple, and CLVD tensor components (Section 3.1 of Vavrycuk, 2015).
  mISO=(1/3)*sum(m);
  mDC=(1/2)*(m(1)-m(3)-abs(m(1)+m(3)-2*m(2)));
  mCLVD=(2/3)*(m(1)+m(3)-2*m(2));
  
  pISO=mISO/sum([abs(mISO),mDC,abs(mCLVD)]);
  pDC=mDC/sum([abs(mISO),mDC,abs(mCLVD)]);
  pCLVD=mCLVD/sum([abs(mISO),mDC,abs(mCLVD)]);
  
  % Calculate Vp/Vs ratio for tensile and shear-tensile sources (Vavrycuk, 2015).
  VpVs=sqrt(((pISO/pCLVD)+1)*(4/3));
  
  % Decompose moment tensor into isotropic, and deviatoric components (Jost & Hermann, 1989; Vavrycuk, 2015).
  Miso=mISO*eye(3);
  Mdev=M-Miso;
  Mdc=mDC*V*[1,0,0;0,0,0;0,0,-1]*V';
  if(pCLVD>0)
      Mclvd=abs(mCLVD)*0.5*V*[2,0,0;0,-1,0;0,0,-1]*V';
  else
      Mclvd=abs(mCLVD)*0.5*V*[1,0,0;0,1,0;0,0,-2]*V';
  end;
  
  % Hudson respresentation variables (Hudson et al., 1989).
  mdev=m-mISO;
  [~,I]=sort(abs(mdev),'descend');
  T_hudson=-mdev(3)/abs(mdev(1));
  k_hudson=mISO/(abs(mISO)+abs(mdev(1)));
  
  % Determine vectors for the P-, B-, and T-axes (Jost & Hermann, 1989).
  [~,I1]=max(mdev);  Taxis=V(:,I1);  Taxis=(Taxis'*(M*Taxis))*Taxis; 
  [~,I3]=min(mdev);  Paxis=V(:,I3);  Paxis=(Paxis'*(M*Paxis))*Paxis; 
  Taxis=Taxis/norm(Taxis);
  Paxis=Paxis/norm(Paxis);
  Baxis=cross(Taxis, Paxis);
  
  % Make sure P- and T-axes are lower hemisphere.
  if(Taxis(3)<0)
      Taxis=-Taxis;
  end;
  if(Paxis(3)<0)
      Paxis=-Paxis;
  end;
  
  % Determine slip vectors (u1,u2) and fault normal vectors (nu1,nu2) (Jost & Hermann, 1989).
  u1 =(Taxis+Paxis)/sqrt(2);
  nu1=(Taxis-Paxis)/sqrt(2);
  u2 =(Taxis-Paxis)/sqrt(2);
  nu2=(Taxis+Paxis)/sqrt(2);
  
  % Make sure fault/auxillary normals are pointing upwards.
  if(nu1(3)>0)
      u1=-u1;
      nu1=-nu1;
  end;
  if(nu2(3)>0)
      u2=-u2;
      nu2=-nu2;
  end;
  
  % Determine fault/auxillary strikes, and dips.
  dip1=acosd(-nu1(3));
  dip2=acosd(-nu2(3));
  strike1=acosd(nu1(2)/sind(dip1));
  strike2=acosd(nu2(2)/sind(dip2));
  
  % Can't determine strike & rake if dip=0, because of div by zero.
  if(sind(dip1)==0)
      fprintf('Can''t determine strike1/rake1 unambiguously since dip1=0\n');
  end;
  if(sind(dip2)==0)
      fprintf('Can''t determine strike2/rake2 unambiguously since dip2=0\n');
  end;
  
  % Put strike values in the proper quadrant.
  Sstr=-nu1(1)/sind(dip1);
  if ( Sstr < 0 ) 
      strike1=360-strike1;
  end
  Sstr=-nu2(1)/sind(dip2);
  if ( Sstr < 0 ) 
      strike2=360-strike2;
  end
  
  % Determine rake values.
  fhsV1=[cosd(strike1), sind(strike1), 0];
  fhsV2=[cosd(strike2), sind(strike2), 0];
  if(u1(3)>=0) % Note the reversal of sign for up-down sense here.
      rake1=-acosd(fhsV1*u1);
  else
      rake1=acosd(fhsV1*u1);
  end
  if(u2(3)>=0)
      rake2=-acosd(fhsV2*u2);
  else
      rake2=acosd(fhsV2*u2);
  end
  
  % P,T,&B plunge angle.
  Pplunge=asind(abs(Paxis(3))/norm(Paxis));
  Tplunge=asind(abs(Taxis(3))/norm(Taxis));
  Bplunge=asind(abs(Baxis(3))/norm(Baxis));
  
  % P-axis azimuth angle, of the two possibilities the positive angle is output.
  Pazimuth=atan2d(Paxis(2),Paxis(1)); 
  if(Pazimuth<0)
      Pazimuth=Pazimuth+360;
  end;
  
  % T-axis azimuth angle, of the two possibilities the positive angle is output.
  Tazimuth=atan2d(Taxis(2),Taxis(1));
  if(Tazimuth<0)
      Tazimuth=Tazimuth+360;
  end;

  % B-axis azimuth angle, of the two possibilities the positive angle is output.
  Bazimuth=atan2d(Baxis(2),Baxis(1));
  if(Bazimuth<0)
      Bazimuth=Bazimuth+360;
  end;
  
  % Define fractions of RF, NR, & SS for DC component (Frohlich, 1992).
  fRF=Taxis(3)^2;
  fNF=Paxis(3)^2;
  fSS=Baxis(3)^2;
  
  % Plot and output information only if flagged to.
  if(plot_flag==1)
      
      % Output information
      fprintf('\n');
      fprintf('Mw: %0.2f\n',Mw);
      fprintf('Mo: %0.2e\n\n',Mo);
      fprintf('pISO:  %4.1f%%\n', pISO*100);
      fprintf('pDC:   %3.1f%%\n', pDC*100);
      fprintf('pCLVD: %3.1f%%\n\n', pCLVD*100);
      fprintf('strike:  %4.1f  %4.1f\n',strike1, strike2);
      fprintf('dip:     %4.1f  %4.1f\n',dip1,dip2);
      fprintf('rake:   %+4.1f %+4.1f\n\n',rake1,rake2);
      fprintf('P,T,&B-azi:    %4.1f %4.1f %4.1f\n',Pazimuth,Tazimuth,Bazimuth);
      fprintf('P,T,&B-plunge:  %4.1f  %4.1f %4.1f\n\n',Pplunge,Tplunge,Bplunge);
      fprintf('pRF:  %4.1f%%\n', fRF*100);
      fprintf('pNF:  %4.1f%%\n', fNF*100);
      fprintf('pSS:  %4.1f%%\n\n', fSS*100);
      
      fprintf('(%3.0f %3.0f %3.0f)  (%3.0f %3.0f %3.0f)  (P: %3.0f %3.0f)  (T: %3.0f %3.0f)\n', strike1, dip1, rake1, strike2, dip2, rake2, Pazimuth, Pplunge, Tazimuth, Tplunge);
      
      % Prep for radiation patterns.
      theta=0:pi/150:2*pi;
      phi=0:pi/150:pi;
      T=[];
      P=T;
      Rp=T;
      Rsv=T;
      Rsh=T;
      
      % Loop through unit sphere.
      for i=1:length(phi);
          for j=1:length(theta)
              
              T=[T,theta(j)];
              P=[P,phi(i)];
              
              % Radiation pattern definitions (Frohlich, 1996; Aki & Richards, 2002).
              % P,Sv,Sh coordinate handedness follows Frohlich (1996).
              unitP=[ sin(phi(i))*cos(theta(j));  sin(phi(i))*sin(theta(j));  cos(phi(i)) ];
              unitSv=[ sin(phi(i)+pi/2)*cos(theta(j));  sin(phi(i)+pi/2)*sin(theta(j));  cos(phi(i)+pi/2) ];
              unitSh=cross(unitSv, unitP);
              
              ampP=unitP'*(M*unitP);
              ampSv=unitSv'*(M*unitP);
              ampSh=unitSh'*(M*unitP);
              
              Rp=[Rp,ampP];
              Rsv=[Rsv,ampSv];
              Rsh=[Rsh,ampSh];
              
              %fprintf('%0.3g %0.3g\n', ampP, ampSh);
          end;
      end;
  
      % Prep Fig1.
      figure(1); clf;
       
      % Plot P-wave radiation pattern.
      subplot(221);
      surf(reshape(abs(Rp).*sin(P).*cos(T),length(theta),[]),reshape(abs(Rp).*sin(P).*sin(T),length(theta),[]),reshape(abs(Rp).*cos(P),length(theta),[]), reshape(Rp,length(theta),[]), 'EdgeColor', 'none' ); hold on;
      scatter3([Paxis(1),-Paxis(1)]*c,[Paxis(2),-Paxis(2)]*c,[Paxis(3),-Paxis(3)]*c, 200,'d','MarkerEdgeColor','k', 'MarkerFaceColor', 'w');
      scatter3([Taxis(1),-Taxis(1)]*c,[Taxis(2),-Taxis(2)]*c,[Taxis(3),-Taxis(3)]*c, 200,'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      scatter3([Baxis(1),-Baxis(1)]*c,[Baxis(2),-Baxis(2)]*c,[Baxis(3),-Baxis(3)]*c, 200,'x','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      title('P-wave radiation pattern'); xlabel('North-South'); ylabel('East-West'); zlabel('Down-Up'); colormap(flipud(parula)); colorbar();
      set(gca, 'Ydir', 'reverse', 'Zdir', 'reverse', 'CLim', clim);
      
      % Plot double couple focal mechanism.
      subplot(222);
      bb([strike1,dip1,rake1],0,0,0,0,'k');
      %bb([Mxx,Myy,Mzz,Mxy,Mxz,Myz],0,0,0,0,'k');
      title('Double couple solution');
      
      % Plot Sv-wave radiation pattern.
      subplot(223);
      surf(reshape(abs(Rsv).*sin(P).*cos(T),length(theta),[]),reshape(abs(Rsv).*sin(P).*sin(T),length(theta),[]),reshape(abs(Rsv).*cos(P),length(theta),[]), reshape(Rsv,length(theta),[]), 'EdgeColor', 'none' ); hold on;
      scatter3([Paxis(1),-Paxis(1)]*c,[Paxis(2),-Paxis(2)]*c,[Paxis(3),-Paxis(3)]*c, 200,'d','MarkerEdgeColor','k', 'MarkerFaceColor', 'w');
      scatter3([Taxis(1),-Taxis(1)]*c,[Taxis(2),-Taxis(2)]*c,[Taxis(3),-Taxis(3)]*c, 200,'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      scatter3([Baxis(1),-Baxis(1)]*c,[Baxis(2),-Baxis(2)]*c,[Baxis(3),-Baxis(3)]*c, 200,'x','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      title('SV-wave radiation pattern'); xlabel('North-South'); ylabel('East-West'); zlabel('Down-Up'); %colorbar();
      set(gca, 'Ydir', 'reverse', 'Zdir', 'reverse', 'CLim', clim);
      
      % Plot Sh-wave radiation pattern.
      subplot(224);
      surf(reshape(abs(Rsh).*sin(P).*cos(T),length(theta),[]),reshape(abs(Rsh).*sin(P).*sin(T),length(theta),[]),reshape(abs(Rsh).*cos(P),length(theta),[]), reshape(Rsh,length(theta),[]), 'EdgeColor', 'none' ); hold on;
      scatter3([Paxis(1),-Paxis(1)]*c,[Paxis(2),-Paxis(2)]*c,[Paxis(3),-Paxis(3)]*c, 200,'d','MarkerEdgeColor','k', 'MarkerFaceColor', 'w');
      scatter3([Taxis(1),-Taxis(1)]*c,[Taxis(2),-Taxis(2)]*c,[Taxis(3),-Taxis(3)]*c, 200,'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      scatter3([Baxis(1),-Baxis(1)]*c,[Baxis(2),-Baxis(2)]*c,[Baxis(3),-Baxis(3)]*c, 200,'x','MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
      title('SH-wave radiation pattern'); xlabel('North-South'); ylabel('East-West'); zlabel('Down-Up'); %colorbar();
      set(gca, 'Ydir', 'reverse', 'Zdir', 'reverse', 'CLim', clim);
      
      % Prep Fig2.
      figure(2); clf;
      
      % Ternary plot: Bottom right thrust, bottom left normal, top strike-slip.
      [hf,vf]=triangle_coords(Pplunge,Tplunge,Bplunge);
      
      Btemp=60*ones(1,101);
      temp=sind(Btemp(1))^2;
      Ptemp=0:(1-temp)/100:(1-temp);
      Ttemp=1-temp-Ptemp;
      %sum([sind(Btemp).^2;Ptemp;Ttemp])
      Ptemp=asind(sqrt(Ptemp));
      Ttemp=asind(sqrt(Ttemp));
      [hSSf,vSSf]=triangle_coords(Ptemp,Ttemp,Btemp);
      
      Ptemp=60*ones(1,101);
      temp=sind(Ptemp(1))^2;
      Btemp=0:(1-temp)/100:(1-temp);
      Ttemp=1-temp-Btemp;
      Btemp=asind(sqrt(Btemp));
      Ttemp=asind(sqrt(Ttemp));
      [hNFf,vNFf]=triangle_coords(Ptemp,Ttemp,Btemp);
      
      Ttemp=50*ones(1,101);
      temp=sind(Ttemp(1))^2;
      Ptemp=0:(1-temp)/100:(1-temp);
      Btemp=1-temp-Ptemp;
      Ptemp=asind(sqrt(Ptemp));
      Btemp=asind(sqrt(Btemp));
      [hRFf,vRFf]=triangle_coords(Ptemp,Ttemp,Btemp);
      
      subplot(221);
      plot(hf,vf, 'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'b'); hold on;
      plot(hSSf,vSSf); % SS bounds.
      plot(hNFf,vNFf); % NF bounds.
      plot(hRFf,vRFf); % RF bounds.
      line([0 -0.5774],[2/3 -1/3],'color','k', 'LineStyle', '--'); % SS-NF boundary.
      line([0 0.5774],[2/3 -1/3],'color','k', 'LineStyle', '--'); % SS-RF boundary.
      line([-0.5774 0.5774],[-1/3 -1/3],'color','k', 'LineStyle', '--'); % NF-RF boundary.
      xlim([-0.65 0.65]); ylim([-0.5 0.8]);
      title('Triangle diagram (Frohlich, 1992)');
      
      
      % Plot: Bottom right is NF, bottom left is SS, and top left is RF.
      subplot(222);
      plot(Pplunge, Tplunge, 'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'b'); hold on;
      line([52 90],[35 35],'color','k', 'LineStyle', '--'); % NF boundary 1.
      line([52 52],[ 0 35],'color','k', 'LineStyle', '--'); % NF boundary 2.
      line([40 52],[20 20],'color','k', 'LineStyle', '--'); % Transtension boundary 1.
      line([40 40],[ 0 20],'color','k', 'LineStyle', '--'); % Transtension boundary 2.
      line([20 40],[20 20],'color','k', 'LineStyle', '--'); % SS boundary 1.
      line([20 20],[20 40],'color','k', 'LineStyle', '--'); % SS boundary 2.
      line([20 20],[40 52],'color','k', 'LineStyle', '--'); % Transpression boundary 1.
      line([ 0 20],[40 40],'color','k', 'LineStyle', '--'); % Transpression boundary 2.
      line([35 35],[52 90],'color','k', 'LineStyle', '--'); % RF boundary 1.
      line([ 0 35],[52 52],'color','k', 'LineStyle', '--'); % RF boundary 2.
      title('Stress Regime Characterization (After Zoback (1992) criteria)'); xlabel('P-axis plunge'); ylabel('T-axis plunge');
      xlim([0 90]); ylim([0 90]);
      
      % Plot CLVD-ISO diagram: middle is double couple, 1st quadrant (top-right) tensile crack, 3rd quadrant compressive crack, other quadrants are withour physical interpretation.
      subplot(223);
      plot(pCLVD, pISO, 'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'b'); hold on;
      line([0 1],[1 0],'color','k');
      line([1 0],[0 -1],'color','k');
      line([0 -1],[-1 0],'color','k');
      line([-1 0],[0 1],'color','k');
      line([0 0],[-1 1],'color','k', 'LineStyle', '--');
      line([-1 1],[0 0],'color','k', 'LineStyle', '--');
      xlim([-1 +1]); ylim([-1 +1]);
      title('CLVD-ISO Plot'); xlabel('CLVD component'); ylabel('ISO component');
      
      % Plot Hudson representation of moment tensor: left is +CLVD, +dipole, and +crack; top is explosion; middle is DC.
      subplot(224);
      plot(T_hudson, k_hudson, 'o','MarkerEdgeColor','k', 'MarkerFaceColor', 'b');
      xlim([-1.5 +1.5]); ylim([-1 +1]);
      title('Hudson representation'); xlabel('-2*epsilon'); ylabel('kappa');
      
  end;


return;
%%%% END OF MAIN FXN.




function [h,v]=triangle_coords(P,T,B)
% Simple function to output triangle diagram coordinates for given P,T,&B plunges.
% Maps points to a ternary diagram based on a azimuthal gnomonic projection, see Frohlich (2001) for details.
% Frohlich, C. (2001). Display and quantitative assessment of distributions of earthquake focal mechanisms. Geophysical Journal International, 144(2), 300-308, doi: 10.1046/j.1365-246x.2001.00341.x.
  
  % Simple projection.
  hs=(sind(T).^2-sind(P).^2)/sqrt(3);
  vs=sind(B).^2-1/3;
  
  % Azimuthal gnomonic projection, from Cliff's Fortran code.
  ct=sind(T);
  cp=sind(P);
  cb=sind(B);
  
  R=sqrt(2)/3;
  s3=sqrt(1/3);
  c3=sqrt(2/3);
  sblat=cb;
  cblat=sqrt(1-cb.*cb);
  plon=(atan2(ct,cp))-pi/4;
  cplon=cos(plon);
  splon=sin(plon);
  
  hg=R.*cblat.*splon./(s3.*sblat+c3.*cblat.*cplon);
  vg=R.*(c3.*sblat-s3.*cblat.*cplon)./(s3.*sblat+c3.*cblat.*cplon);
  
  % Combine the two projections to reduce areal growth and distortion (Frohlich, 2001).
  h=(2/3)*hg+(1/3)*hs;
  v=(2/3)*vg+(1/3)*vs;
  
return;