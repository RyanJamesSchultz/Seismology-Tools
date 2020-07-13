function [dist, azimuth_start, azimuth_end] = Geoid_Distance(lat_start, lon_start, lat_end, lon_end, method)
  % Given a pair of earth coordinates this routine will calculate the
  % surface distance in between them and the azimuth of the great-circle 
  % path at the beginning and ending coordinates.  Lattitude and longitude 
  % are expected as inputs ranging from (-90 to 90) and (-180 to 180).
  %
  % Written by Ryan Schultz.
  
  % Convert coordinates to radians for trig fxns.
  lat_start=lat_start*(pi()/180);
  lat_end=lat_end*(pi()/180);
  lon_start=lon_start*(pi()/180);
  lon_end=lon_end*(pi()/180);
  
  % Azimuth and distance outputs are in degrees.
  if( strcmp(method, 'spherical') )
      [dist, azimuth_start, azimuth_end]=Haversine(lat_start, lon_start, lat_end, lon_end);
      
  elseif( strcmp(method, 'elliptical') )
      [dist, azimuth_start, azimuth_end, status]=Vincenty(lat_start, lon_start, lat_end, lon_end);
      if( status==0 )
          % Antipodal input points cause the code to fail to converge.
          % Dunno how to best handle this atm, just defaulting to shperical if it does happen.
          [dist, azimuth_start, azimuth_end]=Haversine(lat_start, lon_start, lat_end, lon_end);
      end;
      
  else
      fprintf('Enter a proper method: spherical or elliptical\n');
      dist=NaN;
      azimuth_start=NaN;
      azimuth_end=NaN;
      return;
  end;
  
return;
% END OF MAIN FUNCTION.




function [dist, azi_start, azi_end] = Haversine(lat1, lon1, lat2, lon2)
  % Spherical approximation of the earth.
  % Uses Haversine formula to perform spherical trig, ie find great-circle path.
  % Reasonably accurate typically 0.3% error from real earth with up to 0.55%
  % error when crossing the equator.
  
  dLat=lat2-lat1;
  dLon=lon2-lon1;
  
  a = sin(dLat/2.0).*sin(dLat/2.0)+cos(lat1).*cos(lat2).*sin(dLon/2.0).*sin(dLon/2.0);
  c = 2.0*atan2( sqrt(a), sqrt((1.0-a)) );
  
  y = sin(dLon).*cos(lat2);
  x = cos(lat1).*sin(lat2)-sin(lat1).*cos(lat2).*cos(dLon);
  theta = atan2( y, x );
  
  y = sin(-dLon).*cos(lat1);
  x = cos(lat2).*sin(lat1)-sin(lat2).*cos(lat1).*cos(-dLon);

  theta_e = atan2( y, x );
  
  dist=c.*(180/pi()); % Distance in degrees.
  azi_start=rem( (theta*180.0/pi())+360.0, 360.0 );
  azi_end=rem( (theta_e*180.0/pi())+540.0, 360.0 );
  
  
return;



function [dist, azi_start, azi_end, status] = Vincenty(lat1, lon1, lat2, lon2)
  % Uses Vincenty formula to find the distance, coupled with the major and minor 
  % axes measurments from the WGS-84 earth model.  
  % Routine terminates when the error has converged to an accuracy better than a millimetre.
  % Vincenty's method converges slowly or fails for some near anti-podal inputs.
  
  Radius=6371.0;
  a = 6378.137;       b = 6356.752314245;    % WGS-84
 %a = 6378.137;       b = 6356.752314140;    % GRS-80
 %a = 6377.563396;    b = 6356.256910;       % Airy 1830
 %a = 6378.388;       b = 6356.911946;       % International 1924
 %a = 6378.249145;    b = 6356.51486955;     % Clarke model 1880
 %a = 6378.160;       b = 6356.774719;       % GRS-67
  f = (a-b)/a;
  L = lon2-lon1;
  
  U1 = atan( (1.0-f).*tan(lat1) );
  U2 = atan( (1.0-f).*tan(lat2) );
  SinU1 = sin(U1); CosU1 = cos(U1);
  SinU2 = sin(U2); CosU2 = cos(U2);
  lambda = L;
  
  temp=0;
  iter_count=0;
  delta=Inf;
  
  while( (max(delta)>10E-12)&&(iter_count<=10000) )
      SinL = sin(lambda);
      CosL = cos(lambda);
      ss1 = CosU2.*SinL;
      ss2 = CosU1.*SinU2-SinU1.*CosU2.*CosL;
      SinSig = sqrt( ss1.*ss1+ss2.*ss2 );
      CosSig = SinU1.*SinU2+CosU1.*CosU2.*CosL;
      sig = atan2(SinSig, CosSig);
      SinAl = CosU1.*CosU2.*SinL./SinSig;
      CosSqAl = 1.0-SinAl.*SinAl;
      Cos2Sigm = CosSig-2.0.*SinU1.*SinU2./CosSqAl;
      if( max(isnan(Cos2Sigm))||max(isinf(Cos2Sigm)) )
          Cos2Sigm( isnan(Cos2Sigm)|isinf(Cos2Sigm) )=0.0; % On the equator, set to zero.
      end;
      C = (f.*CosSqAl/16.0).*(4.0+f.*(4.0-3.0*CosSqAl));
      delta=temp;
      temp = (1-C).*f.*SinAl.*( sig + C.*SinSig.*( Cos2Sigm + C.*CosSig.*(-1.0 + 2.0*Cos2Sigm.*Cos2Sigm)));
      delta = abs( abs(temp) - abs(delta) );
      lambda = L + temp;
      iter_count=iter_count+1;
  end;
  if(iter_count==10000)
      status=0; % Code failed to converge.
  else
      status=1;
  end;
  
  SinL = sin(lambda);
  CosL = cos(lambda);
  
  uSq = CosSqAl*(a*a-b*b)/(b*b);
  A = 1.0 + (uSq/16384.0).*( 4096.0 + uSq.*( -768.0 + uSq.*( 320.0 - 175.0*uSq )));
  B = (uSq/1024.0).*( 256.0 + uSq.*( -128.0 + uSq.*( 74.0 - 47.0*uSq)));
  ds1 = CosSig.*(-1.0 + 2.0*Cos2Sigm.*Cos2Sigm);
  ds2 = (B/6.0).*Cos2Sigm.*(-3.0 + 4.0*SinSig.*SinSig).*(-3.0 + 4.0*Cos2Sigm.*Cos2Sigm);
  DeltaSig = B.*SinSig.*( Cos2Sigm + 0.25*B.*(  ds1 - ds2  ));
  s = b*A.*(sig-DeltaSig);
  alpha1 = atan2( CosU2.*SinL,  CosU1.*SinU2-SinU1.*CosU2.*CosL );
  alpha2 = atan2( CosU1.*SinL, -SinU1.*CosU2+CosU1.*SinU2.*CosL );
  
  s(isnan(s))=0;
  
  dist=s*180.0/(pi()*Radius);
  azi_start=rem( (alpha1*180.0/pi())+360.0, 360.0);
  azi_end=rem( (alpha2*180.0/pi())+360.0, 360.0);
  
return;
