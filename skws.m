%-------------------------------------------------------------------------%
%               UFSC - Federal University of Santa Catarina               %
%               Graduate Program in Mechanical Engineering                %
%                                                                         %
%     Programmer:                                                         %
%       Davi Klein                                                        %
%                                                                         %
%   Version: 1.0                                              08/09/2022  %
%=========================================================================%
%                          Program Descriprion                            %
%=========================================================================%
%	Function file responsible for the calculation of a skew-symmetric     %
%   matrix of a vector                                                    %
%-------------------------------------------------------------------------%

function sx = skws(S)

sx = [0 -S(3) S(2) ; S(3) 0 -S(1) ; -S(2) S(1) 0 ];

end