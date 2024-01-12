function tono = tonoangle(betas1,betas2)
% output in radians
r           = [zscore(betas1)' zscore(betas2)'];
rUnitCircle = r./sqrt(sum(r.^2,2));
% rad2deg( atan2( sin(pi/4), -cos(pi/4)) ); % for understanding
% interpretations
%                           High                % Low
% tono  = rad2deg( atan2( rUnitCircle(:,1), rUnitCircle(:,2)  ) );
% tono1 = rad2deg( atan( rUnitCircle(:,1)./rUnitCircle(:,2)  ) ); between
% pi/2 and -pi/2
  tono  =   ( atan2( rUnitCircle(:,1), rUnitCircle(:,2)  ) );
% tono  =   ( atan( rUnitCircle(:,1)./rUnitCircle(:,2)  ) );
%                     
%                   High + 90
%                   |
%     180           |
%   Low-   ------------------Low + 0
%                   |
%                   |
%                 High  -  
%                  270
%
%     
%
%