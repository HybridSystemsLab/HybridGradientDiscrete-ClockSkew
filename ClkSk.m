classdef ClkSk2 < HybridSystem
    
    properties(SetAccess = immutable)
        A
        B
        H
        M
        s
        gammac
        gammad
    end

    methods 
        function this = ClkSk(parameters)
            state_dim = 7; % (z(tau_c,q),tauS,tau,thetahat)
            this = this@HybridSystem(state_dim);
            
            this.A = parameters.A;
            this.B = parameters.B;
            this.H = parameters.H;
            this.M = parameters.M;
            this.s = parameters.s;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
        end

        function xdot = flowMap(this, x, t, j)
            z = x(1:2);
            tauS = x(3);
            tau = x(4);
            thetahat = x(5);

            u = 1;
            zdot = this.A*z + this.B*u;
            tauSdot = 1;
            taudot = 1;
            thetahatdot = 0*thetahat;

            %xdot = tauSdot;
            %xdot = [zdot; tauSdot; taudot];
            xdot = [zdot; tauSdot; taudot; thetahatdot; thetahatdot; thetahatdot];
        end

        function xplus = jumpMap(this, x, t, j)
            z = x(1:2);
            tauS = x(3);
            tau = x(4);
            thetahatc = x(5);
            thetahatd = x(6);
            thetahatcd = x(7);

            inD1 = tauS >= this.s;
            inD2 = tau >= 0.5;

            if inD1 % take a sample during flows
                zplus = z;
                tauSplus = 0;
                tauplus = tau;

                psi = tau*z(2);
                y = z(1);
                thetahatplusc = thetahatc + this.s*this.gammac*psi*(y - psi'*thetahatc);
                thetahatplusd = thetahatd;
                thetahatpluscd = thetahatcd + this.s*this.gammac*psi*(y - psi'*thetahatcd);
            elseif inD2  % jump due to the jump set
                tauSplus = tauS;
                tauplus = 0;
                u = 1;
                zplus = this.H*z + this.M*u;
                
                psi = tauplus;
                y = zplus(1);
                thetahatplusc = thetahatc;
                thetahatplusd = thetahatd + this.gammad*psi*(y - psi'*thetahatd)/(1 + this.gammad*norm(psi)^2);
                thetahatpluscd = thetahatcd + this.gammad*psi*(y - psi'*thetahatcd)/(1 + this.gammad*norm(psi)^2);

            else
                zplus = z;
                thetahatplusc = thetahatc;
                thetahatplusd = thetahatd;
                thetahatpluscd = thetahatcd;
                tauSplus = tauS;
                tauplus = tau;
            end

            %xplus = tauSplus;
            %xplus = [zplus; tauSplus; tauplus];
            xplus = [zplus; tauSplus; tauplus; thetahatplusc; thetahatplusd; thetahatpluscd];
        end
        
        function inC = flowSetIndicator(this, x)
            inC = 1;
        end

        function inD = jumpSetIndicator(this, x)
            tauS = x(3);
            tau = x(4);
            
            inD1 = tauS >= this.s;                   % take a sample during flows
            inD2 = tau >= 0.5;    % jump due to jump set
            inD = inD1 || inD2;
        end
    end
end
