classdef rateandstate2D < unicycle.ode.evolution
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=rateandstate2D(varargin)
            obj=obj@unicycle.ode.evolution(varargin{:});
            obj.flt.dgf=4;
        end % constructor
        
        function yp=ode(o,t,y)
            % ODE is a method of class RATESTRENGTHENINGPOWER used to
            % solve the quasi-static equations for fault slip evolution
            % when slip is governed by rate-and-state friction:
            %
            %   [t,y]=ode.ode45(@evl.ratestate,[0 200],y0,options);
            %
            % where evl is an instance of class evolution.
            %
            % the rate-and-state friction consistutive equations are:
            %
            %       ds
            %   V = -- = 2*Vo*sinh(tau/(a*sigma))
            %       dt
            %
            % and
            %
            %   d theta
            %   ------- = (Vo*exp(-theta)-V)/L
            %      dt
            %
            % where a, b, Vo, and L are friction properties, and theta
            % is the state variable that describes the healing and weakening
            % of the fault contact.
            %
            %        /        s          \
            %        |       tau         |
            %  y =   | log(theta Vo / L) |
            %        \   log( V / Vo )   /
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            %
            %                     F A U L T S
            %
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            % initialize yp
            yp=zeros(size(y));
            
            % Fault tractions
            tau = y(2:o.flt.dgf:(o.flt.N*o.flt.dgf));
            
            % State variable
            th = y(3:o.flt.dgf:(o.flt.N*o.flt.dgf));
            theta = o.flt.l.*exp(th)./o.flt.Vo;
            
            % Slip velocity
            V    = o.flt.Vo.*exp(y(4:o.flt.dgf:(o.flt.N*o.flt.dgf)));
            
            % pin patches
            if 0~=numel(o.flt.pinnedPosition)
                tau(o.flt.pinnedPosition)=0;
                V(o.flt.pinnedPosition)=0;
            end
            
            % velocity
            yp(1:o.flt.dgf:(o.flt.N*o.flt.dgf))=V;
            
            % Rate of state (rate of log(theta/theta0))
            dth=(o.flt.Vo.*exp(-th)-V)./o.flt.l;
            yp(3:o.flt.dgf:(o.flt.N*o.flt.dgf))=dth;
            
            % Acceleration (rate of log(V/Vo)) - STRIKE SLIP ONLY
            kv = o.KK{1,1}*(V-o.flt.Vpl);
            
            
            % Radiation damping
            damping=0.2.*o.flt.earthModel.G./o.flt.Vs/2;
            
            % Acceleraton
            sigmat=o.flt.sigma;
            yp(4:o.flt.dgf:(o.flt.N*o.flt.dgf)) = (o.flt.mu0./o.flt.a).*(kv./tau)-(o.flt.b./o.flt.a).*dth;            
            yp(4:o.flt.dgf:(o.flt.N*o.flt.dgf)) = yp(4:o.flt.dgf:(o.flt.N*o.flt.dgf))./(o.flt.a.*sigmat + damping.*V);
            
            % Shear stress rate on fault - adding damping
            yp(2:o.flt.dgf:(o.flt.N*o.flt.dgf)) = kv - damping.*V;
            
        end % end ratestate
        
    end % end methods
    
end
