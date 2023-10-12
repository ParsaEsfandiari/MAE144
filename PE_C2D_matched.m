function [Dz]=PE_C2D_matched(Ds,h,causal,omegabar)
        % function PE_C2D_matched(Ds,h,omegabar,causal)
        % Convert Ds(s) to Dz(z) using matched method.  If omegabar is not
        % specified, omegabar is set to be 0. If causal==1, the output is
        % going to be causal and if not specified, the output is going to
        % be non casual.
        % TEST:
        %   h=0.01; Ds=RR_tf([1,1],[1,11,1]); PE_C2D_matched(Ds,h,1)
        %   disp('Corresponding Matlab solution:')
        %   s=tf('s'); Dc=(s+1)/((s+1)*(s+10));c2d(Dc,0.01,'matched')
        % TEST for symbolic solution:
        % syms z1 p1; h=0.01; Ds=RR_tf([1,z1],[1,p1,1]); PE_C2D_matched(Ds,h,1)
        % Parsa Esfandiari, hw1, https://github.com/ParsaEsfandiari/
            if nargin==2, omegabar=0; causal=0; 
            elseif nargin==3, omegabar=0; end
            
            m=Ds.num.n; n=Ds.den.n; b=RR_poly(0); a=b;
            N=roots(Ds.num.poly); M=roots(Ds.den.poly);
            N=exp(N*h);M=exp(M*h);
            
            SGainA=RR_poly(0); SGainB=RR_poly(0);
            if omegabar==0
                SGainA=Ds.den.poly(n+1);
                SGainB=Ds.num.poly(m+1);
            else
            for j=1:n+1; SGainA=SGainA+Ds.den.poly(j)*(omegabar^(n-j+1));    end
            for j=1:m+1; SGainB=SGainB+Ds.num.poly(j)*(omegabar^(m-j+1));    end
            end
            
            if causal~=1 
                for j=1:1:n-m,
                N=[N,-1];
                end
                b=RR_poly([N],1);
            else
                for j=1:1:n-m-1,
                    N=[N,-1];
                end
                b=RR_poly([N],1);
            end
            a=RR_poly([M],1);
            
            A=RR_tf(a); B=RR_tf(b);
            ZGainA=RR_poly(0); ZGainB=RR_poly(0);
            if omegabar==0
                for j=1:(A.num.n)+1; ZGainA=ZGainA+A.num.poly(j);    end
                for j=1:(B.num.n)+1; ZGainB=ZGainB+B.num.poly(j);    end
            else
                for j=1:(A.num.n)+1; ZGainA=ZGainA+A.num.poly(j)*(exp(omegabar*h));    end
                for j=1:(B.num.n)+1; ZGainB=ZGainB+B.num.poly(j)*(exp(omegabar*h));    end
            end
            K=RR_poly(0);
            k=(ZGainA*SGainB)/(ZGainB*SGainA);
            Dz=RR_tf(k*b,a);
        end % function PE_C2D_matched
