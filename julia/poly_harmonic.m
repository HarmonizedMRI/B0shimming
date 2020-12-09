function out = poly_harmonic(n, m, pos_x, pos_y, pos_z)
% returns harmonic field for the required solid harmonic addressed by
% n, m. Positive m values correspond to cosine component and negative
% to sine
% pos_... can be both value and vector/matrix

r2=pos_x.^2+pos_y.^2; %pos_z.^2 is out!

if n==0,
	out=1;
else if n==1,
	if m==0,
		out=pos_z;
	else
		if m==1,
			out=-pos_x;
		else
			out=-pos_y; %m==-1
		end
	end
else if n==2,
	if m==0,
		out=pos_z.^2-0.5*r2;
	else if abs(m)==1,
		if m==1,
			out=-pos_z.*pos_x;
		else
                        out=-pos_z.*pos_y; %m==-1
                    end
                else
                    if m==2,
                        out=pos_x.^2-pos_y.^2;
                    else
                        out=2*pos_x.*pos_y; %m==-2
                    end
                end
            end
        else if n==3,
                if m==0,
                    out=pos_z.*(pos_z.^2-2/3*r2);
                else if abs(m)==1,
                        if m==1,
                            out=-pos_x.*(pos_z.^2-0.25*r2);
                        else
                            out=-pos_y.*(pos_z.^2-0.25*r2);%m==-1
                        end
                    else if abs(m)==2,
                            if m==2,
                                out=pos_z.*(pos_x.^2-pos_y.^2);
                            else
                                out=2*pos_z.*pos_x.*pos_y;%m==-2
                            end;
                        else
                            if m==3,
                                out=-pos_x.*(pos_x.^2-3*pos_y.^2);
                            else
                                out=-pos_y.*(3*pos_x.^2-pos_y.^2);%m==-
3
                            end
                        end
                    end
                end
            else if n==4,
                    if m==0,
                        out=pos_z.^4-3*pos_z.^2.*r2+3/8*r2.^2;
                    else if abs(m)==1,
                            if m==1,
                                out=-pos_z.*pos_x.*(pos_z.^2-0.75*r2);
                            else
                                out=-pos_z.*pos_y.*(pos_z.^2-0.75*r2);%m==-1

                            end
                        else if abs(m)==2,
                                if m==2,
                                    out=(pos_x.^2-pos_y.^2).*(pos_z.^2-1/6*r2);

                                else
                                    out=2*pos_x.*pos_y.*(pos_z.^2-1/6*r2);%m==-2

                                end
                            else if abs(m)==3,
                                    if m==3,
                                        out=-pos_z.*pos_x.*(pos_x.^2-3*pos_y.^2);

                                    else
                                        out=-pos_z.*pos_y.*(3*pos_x.^2-pos_y.^2);%m==-3

                                    end
                                else
                                    if m==4,
                                        out=pos_x.^4-6*pos_x.^2.*pos_y.^2+pos_y.^4;

                                    else
                                        out=4*pos_x.*pos_y.*(pos_x.^2-pos_y.^2);%m==-4

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
