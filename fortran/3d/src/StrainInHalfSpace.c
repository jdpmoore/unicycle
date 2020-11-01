// inialize the libraries
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <float.h>

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Calculate main dislocation contribution to strains and stresses
void       TDstressFS_inStrnHS(double StsMS[6], double StrMS[6],  double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double mu, double lambda);
//Calculate harmonic function contribution to strains and stresses
void TDstress_HarFunc_inStrnHS(double StsFSC[6],double StrFSC[6], double X, double Y, double Z, double P1[3] ,double P2[3], double P3[3], double Ss, double Ds, double Ts, double mu, double lambda);

void       CoordTrans_inStrnHS(double newVal[3], double x_shift, double y_shift, double z_shift, double RotMat[3][3]); 

void      TriModeFind_inStrnHS(int TrimMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b);

void         TDSetupS_inStrnHS(double x,double y,double z,double alpha,double bx,double by,double bz,double nu, double TriVertex[3],double SideVec[3],double e[6]);

void     AngDisStrain_inStrnHS(double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double e[6]);

void        TensTrans_inStrnHS(double e_in[6], double e_out[6], double B[3][3]);

void    AngSetupFSC_S_inStrnHS(double Stress1[6],double Strain1[6], double X,double Y,double Z,double bX,double bY,double bZ,double Pt1[3], double Pt2[3], double mu,double lambda);

void  AngDisStrainFSC_inStrnHS(double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a, double Strain[6]);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void straininhalfspace_(double Stress[6], double Strain[6], double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double mu, double lambda)
/*
this function is translated by Olaf Zielke from Matlab to C
% TDstressHS 
% Calculates stresses and strains associated with a triangular dislocation 
% in an elastic half-space.
%
% TD: Triangular Dislocation
% EFCS: Earth-Fixed Coordinate System
% TDCS: Triangular Dislocation Coordinate System
% ADCS: Angular Dislocation Coordinate System
% 
% INPUTS
% X, Y and Z: 
% Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
% must have the same size.
%
% P1,P2 and P3:
% Coordinates of TD vertices in EFCS.
% 
% Ss, Ds and Ts:
% TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
%
% mu and lambda:
% Lame constants.
%
% OUTPUTS
% Stress:
% Calculated stress tensor components in EFCS. The six columns of Stress 
% are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components 
% have the same unit as Lame constants.
%
% Strain:
% Calculated strain tensor components in EFCS. The six columns of Strain 
% are Exx, Eyy, Ezz, Exy, Exz and Eyz, respectively. The strain components 
% are dimensionless.
% 
% 
% Example: Calculate and plot the first component of stress tensor on a  
% regular grid.
% 
% [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,-5);
% [Stress,Strain] = TDstressHS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 -2],
% -1,2,3,.33e11,.33e11);
% h = surf(X,Y,reshape(Stress(:,1),size(X)),'edgecolor','none');
% view(2)
% axis equal
% axis tight
% set(gcf,'renderer','painters')

% Reference journal article: 
% Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
% artefact-free solution. 
% Submitted to Geophysical Journal International 

% Copyright (c) 2014 Mehdi Nikkhoo
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

% I appreciate any comments or bug reports.

% Mehdi Nikkhoo
% created: 2013.1.28
% Last modified: 2014.7.30
% 
% VolcanoTectonics Research Group
% Section 2.1, Physics of Earthquakes and Volcanoes
% Department 2, Physics of the Earth
% Helmholtz Centre Potsdam
% German Research Centre for Geosciences (GFZ)
% 
% email: 
% mehdi.nikkhoo@gfz-potsdam.de 
% mehdi.nikkhoo@gmail.com
*/
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
{
    int    i;
    double StsMS[6],                StrMS[6];
    double StsFSC[6],               StrFSC[6];
    double StsIS[6],                StrIS[6];
    double P1mirr[3],               P2mirr[3],                  P3mirr[3];
    
    P1mirr[0] = P1[0];              P2mirr[0] = P2[0];          P3mirr[0] = P3[0];
    P1mirr[1] = P1[1];              P2mirr[1] = P2[1];          P3mirr[1] = P3[1];
    P1mirr[2] = -1.0*P1[2];         P2mirr[2] = -1.0*P2[2];     P3mirr[2] = -1.0*P3[2];
    
    for (i = 0; i < 6; i++)
    {   StsMS[i]   = 0.0;           StrMS[i]  = 0.0;     
        StsFSC[i]  = 0.0;           StrFSC[i] = 0.0;
        StsIS[i]   = 0.0;           StrIS[i]  = 0.0;
    }
    // Calculate main dislocation contribution to strains and stresses
    TDstressFS_inStrnHS(       StsMS, StrMS,  X, Y, Z, P1,     P2,     P3, Ss, Ds, Ts, mu, lambda);

    // Calculate harmonic function contribution to strains and stresses
    TDstress_HarFunc_inStrnHS(StsFSC, StrFSC, X, Y, Z, P1,     P2,     P3, Ss, Ds, Ts, mu, lambda);

    // Calculate image dislocation contribution to strains and stresses
    TDstressFS_inStrnHS(       StsIS, StrIS,  X, Y, Z, P1mirr, P2mirr, P3mirr, Ss, Ds, Ts, mu, lambda);

    if ((P1mirr[2] == 0.0) && (P2mirr[2] == 0.0) && (P3mirr[2] == 0.0))
    {   StsIS[4] = -1.0*StsIS[4];           StsIS[5] = -1.0*StsIS[5];
        StrIS[4] = -1.0*StrIS[4];           StrIS[5] = -1.0*StrIS[5];
    }
    // Calculate the complete stress and strain tensor components in EFCS
    for (i = 0; i < 6; i++)
    {   Stress[i] = StsMS[i] +StsIS[i] +StsFSC[i];
        Strain[i] = StrMS[i] +StrIS[i] +StrFSC[i];
    } 
    
    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void TDstressFS_inStrnHS(double StsMS[6], double StrMS[6],  double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double by, double bz, double bx, double mu, double lambda)
{
    //TDstressFS : Calculates stresses and strains associated with a triangular dislocation in an elastic full-space.
    //bx = Ts; % Tensile-slip;      by = Ss; % Strike-slip;         bz = Ds; % Dip-slip
    int     casepLog,               casenLog;
    double   A,      B,      C,  x,  y,      z,      nu,     Tempdouble;
    
    int      TriMode[1];      
    double   tempVect1[3],   tempVect2[3];
    double   p1[3],          p2[3],			 p3[3];              
    double   e12[3],         e13[3],         e23[3];                 
    double   e_vals1[6],     e_vals2[6],     e_vals3[6],     e_comb[6];
    double   Amat[3][3];
    //--------------------------------------------------------------------------------------------------------------------
    nu = lambda/(lambda+mu) /2.0; // Poisson's ratio
    p1[0]        = 0.0;                     p1[1]        = 0.0;                     p1[2]        = 0.0;
    p2[0]        = 0.0;                     p2[1]        = 0.0;                     p2[2]        = 0.0;
    p3[0]        = 0.0;                     p3[1]        = 0.0;                     p3[2]        = 0.0;
    //--------------------------------//--------------------------------
    //--------------------------------//--------------------------------
    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward, respectively.
    double eY[3],					eZ[3];
    double Vnorm[3],				Vstrike[3],				Vdip[3];
    
    eY[0]        = 0.0;                     eY[1]        = 1.0;                     eY[2]        = 0.0;
    eZ[0]        = 0.0;                     eZ[1]        = 0.0;                     eZ[2]        = 1.0;    
    //--------------------------------
	tempVect1[0] = P2[0] -P1[0];            tempVect1[1] = P2[1] -P1[1];            tempVect1[2] = P2[2] -P1[2];
    tempVect2[0] = P3[0] -P1[0];            tempVect2[1] = P3[1] -P1[1];            tempVect2[2] = P3[2] -P1[2];
      
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
    Tempdouble   = pow (  ( pow(Vnorm[0],2.0) + pow(Vnorm[1],2.0) +pow(Vnorm[2],2.0) ),0.5);
    //--------------------------------
    Vnorm[0]     = Vnorm[0]/Tempdouble;     Vnorm[1]     = Vnorm[1]/Tempdouble;     Vnorm[2]     = Vnorm[2]/Tempdouble;
	//--------------------------------
	Vstrike[0]   = eZ[1]*Vnorm[2] - eZ[2]*Vnorm[1];
    Vstrike[1]   = eZ[2]*Vnorm[0] - eZ[0]*Vnorm[2];
    Vstrike[2]   = eZ[0]*Vnorm[1] - eZ[1]*Vnorm[0];
    // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction)
    Tempdouble   = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
	if (Tempdouble == 0.0)
    {   Vstrike[0] = eY[0]*Vnorm[2];        Vstrike[1] = eY[1]*Vnorm[2];            Vstrike[2]   = eY[2]*Vnorm[2];        
        // For horizontal elements in case of half-space calculation!!! => Correct the strike vector of image dislocation only
        if (P1[2] > 0.0)
        {   Vstrike[0] = -1.0*Vstrike[0];   Vstrike[1] = -1.0*Vstrike[1];           Vstrike[2]   = -1.0*Vstrike[2];
    }   }
	Tempdouble = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    //--------------------------------
    Vstrike[0]    = Vstrike[0]/Tempdouble;  Vstrike[1] = Vstrike[1]/Tempdouble;     Vstrike[2]   = Vstrike[2]/Tempdouble;
    //--------------------------------
    Vdip[0]       = Vnorm[1]*Vstrike[2] - Vnorm[2]*Vstrike[1];
    Vdip[1]       = Vnorm[2]*Vstrike[0] - Vnorm[0]*Vstrike[2];
    Vdip[2]       = Vnorm[0]*Vstrike[1] - Vnorm[1]*Vstrike[0];
    Tempdouble    = pow (  ( pow(Vdip[0],2.0) + pow(Vdip[1],2.0) +pow(Vdip[2],2.0) ),0.5);
    //--------------------------------
    Vdip[0]       = Vdip[0]/Tempdouble;     Vdip[1] = Vdip[1]/Tempdouble;           Vdip[2]      = Vdip[2]/Tempdouble;
    //--------------------------------//--------------------------------
    //--------------------------------//--------------------------------


    // Transform coordinates and slip vector components from EFCS into TDCS
    Amat[0][0]    = Vnorm[0];                   Amat[0][1]    = Vnorm[1];                   Amat[0][2]    = Vnorm[2];
    Amat[1][0]    = Vstrike[0];                 Amat[1][1]    = Vstrike[1];                 Amat[1][2]    = Vstrike[2];
    Amat[2][0]    = Vdip[0];                    Amat[2][1]    = Vdip[1];                    Amat[2][2]    = Vdip[2];    

    /////////////////////////////////////////////////////////////////////////////////////////////    
    CoordTrans_inStrnHS(tempVect1, (X-P2[0]),  (Y-P2[1]), (Z-P2[2]), Amat); 
    x          = tempVect1[0];                  y          = tempVect1[1];                  z          = tempVect1[2];
    CoordTrans_inStrnHS(tempVect1, (P1[0]-P2[0]),  (P1[1]-P2[1]), (P1[2]-P2[2]), Amat); 
    p1[0]      = tempVect1[0];                   p1[1]      = tempVect1[1];                  p1[2]      = tempVect1[2];
    CoordTrans_inStrnHS(tempVect2, (P3[0]-P2[0]),  (P3[1]-P2[1]), (P3[2]-P2[2]), Amat); 
    p3[0]      = tempVect2[0];                   p3[1]      = tempVect2[1];                  p3[2]      = tempVect2[2];
    // Calculate the unit vectors along TD sides in TDCS
    Tempdouble  = pow (  ( pow((p2[0]-p1[0]),2.0) + pow((p2[1]-p1[1]),2.0) +pow((p2[2]-p1[2]),2.0 )),0.5);
    e12[0]     = (p2[0]-p1[0])/Tempdouble;       e12[1]      = (p2[1]-p1[1])/Tempdouble;       e12[2]     = (p2[2]-p1[2])/Tempdouble;     
    Tempdouble  = pow (  ( pow((p3[0]-p1[0]),2.0) + pow((p3[1]-p1[1]),2.0) +pow((p3[2]-p1[2]),2.0 )),0.5);
    e13[0]     = (p3[0]-p1[0])/Tempdouble;       e13[1]      = (p3[1]-p1[1])/Tempdouble;       e13[2]     = (p3[2]-p1[2])/Tempdouble;
    Tempdouble  = pow (  ( pow((p3[0]-p2[0]),2.0) + pow((p3[1]-p2[1]),2.0) +pow((p3[2]-p2[2]),2.0 )),0.5);
    e23[0]     = (p3[0]-p2[0])/Tempdouble;       e23[1]      = (p3[1]-p2[1])/Tempdouble;       e23[2]     = (p3[2]-p2[2])/Tempdouble; 
    // Calculate the TD angles
    Tempdouble  = e12[0]*e13[0] +e12[1]*e13[1] +e12[2]*e13[2];
    A = acos(Tempdouble) ;
    Tempdouble  = -1.0*e12[0]*e23[0] + -1.0*e12[1]*e23[1] + -1.0*e12[2]*e23[2];
    B = acos(Tempdouble) ;
    Tempdouble  = e23[0]*e13[0] +e23[1]*e13[1] +e23[2]*e13[2];
    C = acos(Tempdouble) ;
     
    // Determine the best arteact-free configuration for each calculation point
    ///////////////////////////////////////////////////////////////////////////////////////////// 
    TriModeFind_inStrnHS(TriMode,y,z,x,p1[1],p1[2], p2[1], p2[2], p3[1], p3[2]);
    ///////////////////////////////////////////////////////////////////////////////////////////// 
    if (TriMode[0] == 1)       {       casepLog = 1;   casenLog = 0;          } 
    if (TriMode[0] ==-1)       {       casepLog = 0;   casenLog = 1;          } 
    if (TriMode[0] == 0)       {       casepLog = 0;   casenLog = 0;          } 


    
    if (casepLog == 1) // Configuration I
    {   // Calculate first angular dislocation contribution
        tempVect1[0] = -1.0*e13[0];         tempVect1[1] = -1.0*e13[1];         tempVect1[2] = -1.0*e13[2];              
        TDSetupS_inStrnHS(x, y, z, A, bx, by, bz, nu, p1, tempVect1, e_vals1);
        // Calculate second angular dislocation contribution
        TDSetupS_inStrnHS(x, y, z, B, bx, by, bz, nu, p2, e12,       e_vals2);
        // Calculate third angular dislocation contribution
        TDSetupS_inStrnHS(x, y, z, C, bx, by, bz, nu, p3, e23,       e_vals3);  
    }
    if (casenLog == 1) // Configuration II
    {  // Calculate first angular dislocation contribution             
        TDSetupS_inStrnHS(x, y, z, A, bx, by, bz, nu, p1, e13,       e_vals1);
        // Calculate second angular dislocation contribution
        tempVect1[0] = -1.0*e12[0];         tempVect1[1] = -1.0*e12[1];         tempVect1[2] = -1.0*e12[2];      
        TDSetupS_inStrnHS(x, y, z, B, bx, by, bz, nu, p2, tempVect1, e_vals2);
        // Calculate third angular dislocation contribution
        tempVect1[0] = -1.0*e23[0];         tempVect1[1] = -1.0*e23[1];         tempVect1[2] = -1.0*e23[2]; 
        TDSetupS_inStrnHS(x, y, z, C, bx, by, bz, nu, p3, tempVect1, e_vals3);     
    }
    if ((casenLog == 1) || (casepLog == 1))
    {
        e_comb[0]       = e_vals1[0]+e_vals2[0]+e_vals3[0]; // exx
        e_comb[1]       = e_vals1[1]+e_vals2[1]+e_vals3[1]; // exy
        e_comb[2]       = e_vals1[2]+e_vals2[2]+e_vals3[2]; // exz
        e_comb[3]       = e_vals1[3]+e_vals2[3]+e_vals3[3]; // eyy
        e_comb[4]       = e_vals1[4]+e_vals2[4]+e_vals3[4]; // eyz
        e_comb[5]       = e_vals1[5]+e_vals2[5]+e_vals3[5]; // ezz
    }
    else
    {   
        e_comb[0]       = NAN; // exx => supposed to be "NaN" => have to check
        e_comb[1]       = NAN; // exy
        e_comb[2]       = NAN; // exz
        e_comb[3]       = NAN; // eyy
        e_comb[4]       = NAN; // eyz
        e_comb[5]       = NAN; // ezz
    }
    Amat[0][0] = Vnorm[0];         Amat[0][1] = Vstrike[0];       Amat[0][2] = Vdip[0];
    Amat[1][0] = Vnorm[1];         Amat[1][1] = Vstrike[1];       Amat[1][2] = Vdip[1];
    Amat[2][0] = Vnorm[2];         Amat[2][1] = Vstrike[2];       Amat[2][2] = Vdip[2];
    // Transform the strain tensor components from TDCS into EFCS
    TensTrans_inStrnHS(e_comb,StrMS, Amat);

    // Calculate the stress tensor components in EFCS
    StsMS[0] = 2.0*mu*StrMS[0]+lambda*(StrMS[0]+StrMS[3]+StrMS[5]); // sxx
    StsMS[3] = 2.0*mu*StrMS[3]+lambda*(StrMS[0]+StrMS[3]+StrMS[5]); // syy
    StsMS[5] = 2.0*mu*StrMS[5]+lambda*(StrMS[0]+StrMS[3]+StrMS[5]); // szz
    StsMS[1] = 2.0*mu*StrMS[1]; // sxy
    StsMS[2] = 2.0*mu*StrMS[2]; // sxz
    StsMS[4] = 2.0*mu*StrMS[4]; // syz
    
    return;
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
void TDstress_HarFunc_inStrnHS(double StsFSC[6], double StrFSC[6], double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double by, double bz, double bx, double mu, double lambda)
{
    // TDstress_HarFunc calculates the harmonic function contribution to the strains and stresses associated with a triangular dislocation in a 
    // half-space. The function cancels the surface normal tractions induced by the main and image dislocations.

    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward, respectively.
    int   i;          
    double tempVect1[3],     tempVect2[3];
    double Stress1[6],       Stress2[6],         Stress3[6];
    double Strain1[6],       Strain2[6],         Strain3[6];
    
    double Amat[3][3];
	//--------------------------------//--------------------------------
    //--------------------------------//--------------------------------
    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward, respectively.
    double Tempdouble,				eY[3],					eZ[3];
    double Vnorm[3],				Vstrike[3],				Vdip[3];
    
    eY[0]        = 0.0;                     eY[1]        = 1.0;                     eY[2]        = 0.0;
    eZ[0]        = 0.0;                     eZ[1]        = 0.0;                     eZ[2]        = 1.0;    
    //--------------------------------
    tempVect1[0] = P2[0] -P1[0];            tempVect1[1] = P2[1] -P1[1];            tempVect1[2] = P2[2] -P1[2];
    tempVect2[0] = P3[0] -P1[0];            tempVect2[1] = P3[1] -P1[1];            tempVect2[2] = P3[2] -P1[2];
      
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
    Tempdouble   = pow (  ( pow(Vnorm[0],2.0) + pow(Vnorm[1],2.0) +pow(Vnorm[2],2.0) ),0.5);
    //--------------------------------
    Vnorm[0]     = Vnorm[0]/Tempdouble;     Vnorm[1]     = Vnorm[1]/Tempdouble;     Vnorm[2]     = Vnorm[2]/Tempdouble;
	//--------------------------------
	Vstrike[0]   = eZ[1]*Vnorm[2] - eZ[2]*Vnorm[1];
    Vstrike[1]   = eZ[2]*Vnorm[0] - eZ[0]*Vnorm[2];
    Vstrike[2]   = eZ[0]*Vnorm[1] - eZ[1]*Vnorm[0];
    // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction)
    Tempdouble   = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
	if (Tempdouble == 0.0)
    {   Vstrike[0] = eY[0]*Vnorm[2];        Vstrike[1] = eY[1]*Vnorm[2];            Vstrike[2]   = eY[2]*Vnorm[2];        
    }
	Tempdouble = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    //--------------------------------
    Vstrike[0]    = Vstrike[0]/Tempdouble;  Vstrike[1] = Vstrike[1]/Tempdouble;     Vstrike[2]   = Vstrike[2]/Tempdouble;
    //--------------------------------
    Vdip[0]       = Vnorm[1]*Vstrike[2] - Vnorm[2]*Vstrike[1];
    Vdip[1]       = Vnorm[2]*Vstrike[0] - Vnorm[0]*Vstrike[2];
    Vdip[2]       = Vnorm[0]*Vstrike[1] - Vnorm[1]*Vstrike[0];
    Tempdouble    = pow (  ( pow(Vdip[0],2.0) + pow(Vdip[1],2.0) +pow(Vdip[2],2.0) ),0.5);
    //--------------------------------
    Vdip[0]       = Vdip[0]/Tempdouble;     Vdip[1] = Vdip[1]/Tempdouble;           Vdip[2]      = Vdip[2]/Tempdouble;
    //--------------------------------//--------------------------------
    //--------------------------------//--------------------------------
    
    // Transform slip vector components from TDCS into EFCS
    Amat[0][0] = Vnorm[0];          Amat[0][1] = Vstrike[0];            Amat[0][2] = Vdip[0];
    Amat[1][0] = Vnorm[1];          Amat[1][1] = Vstrike[1];            Amat[1][2] = Vdip[1];
    Amat[2][0] = Vnorm[2];          Amat[2][1] = Vstrike[2];            Amat[2][2] = Vdip[2];

    CoordTrans_inStrnHS(tempVect1, bx,  by, bz, Amat); 
    // Calculate contribution of angular dislocation pair on each TD side 
    AngSetupFSC_S_inStrnHS(Stress1,Strain1, X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P1,P2,mu,lambda); // P1P2   
    AngSetupFSC_S_inStrnHS(Stress2,Strain2, X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P2,P3,mu,lambda); // P2P3
    AngSetupFSC_S_inStrnHS(Stress3,Strain3, X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P3,P1,mu,lambda); // P3P1
    //Calculate total harmonic function contribution to strains and stresses
    for (i = 0; i < 6; i++)
    {   StsFSC[i] = Stress1[i] + Stress2[i] + Stress3[i];
        StrFSC[i] = Strain1[i] + Strain2[i] + Strain3[i];
    }
    
    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void  CoordTrans_inStrnHS(double newVal[3], double x1, double x2, double x3, double A[3][3]) 
{
    // CoordTrans_inStrnHS transforms the coordinates of the vectors, from x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
    // transformation matrix, whose columns e1,e2 and e3 are the unit base vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
    // in X1X2X3. The transpose of A (i.e., A') will transform the coordinates  from X1X2X3 into x1x2x3.
    
    newVal[0] = A[0][0]*x1 + A[0][1]*x2 + A[0][2]*x3;     newVal[1] = A[1][0]*x1 + A[1][1]*x2 + A[1][2]*x3;       newVal[2] = A[2][0]*x1 + A[2][1]*x2 + A[2][2]*x3;

    return;
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void      TriModeFind_inStrnHS(int TriMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b)
{
    // trimodefinder calculates the normalized barycentric coordinates of  the points with respect to the TD vertices and specifies the appropriate
    // artefact-free configuration of the angular dislocations for the  calculations. The input matrices x, y and z share the same size and
    // correspond to the y, z and x coordinates in the TDCS, respectively. p1, p2 and p3 are two-component matrices representing the y and z coordinates
    // of the TD vertices in the TDCS, respectively. The components of the output (trimode) corresponding to each calculation 
    // points, are 1 for the first configuration, -1 for the second configuration and 0 for the calculation point that lie on the TD sides.
    double a,            b,          c;

    a = ((p2_b-p3_b)*(x-p3_a) +(p3_a-p2_a)*(y-p3_b)) /  ((p2_b-p3_b)*(p1_a-p3_a) +(p3_a-p2_a)*(p1_b-p3_b));
    b = ((p3_b-p1_b)*(x-p3_a) +(p1_a-p3_a)*(y-p3_b)) /  ((p2_b-p3_b)*(p1_a-p3_a) +(p3_a-p2_a)*(p1_b-p3_b));
    c = 1.0 -a -b;

    TriMode[0] = 1;
    if       ((a <= 0.0) && (b > c) && (c > a))             {   TriMode[0] = -1;                  }
    else if  ((b <= 0.0) && (c > a) && (a > b))             {   TriMode[0] = -1;                  }
    else if  ((c <= 0.0) && (a > b) && (b > c))             {   TriMode[0] = -1;                  }

    else if  ((a == 0.0) && (b >= 0.0) && (c >= 0.0))       {   TriMode[0] = 0;                   }
    else if  ((a >= 0.0) && (b == 0.0) && (c >= 0.0))       {   TriMode[0] = 0;                   }
    else if  ((a >= 0.0) && (b >= 0.0) && (c == 0.0))       {   TriMode[0] = 0;                   }
    if       ((TriMode[0] == 0) && (z != 0.0))              {   TriMode[0] = 1;                   } 

    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void TDSetupS_inStrnHS(double x,double y,double z,double alpha,double bx,double by,double bz,double nu, double TriVertex[3],double SideVec[3],double e_out[6])
// TDSetupS transforms coordinates of the calculation points as well as slip vector components from ADCS into TDCS. It then calculates the strains in ADCS and transforms them into TDCS.
{
    double A[2][2];
    double B[3][3];
    double r1[2];
    double r2[2];
    double y1;
    double z1;
    double by1;
    double bz1;
    double tempVect1[2];
    double e[6];
    //-----------------------------------
    // Transformation matrix
    A[0][0] = SideVec[2];                   A[0][1]      = -1.0*SideVec[1];
    A[1][0] = SideVec[1];                   A[1][1]      =      SideVec[2]; 
    // Transform coordinates of the calculation points from TDCS into ADCS
    tempVect1[0] = y -TriVertex[1];         tempVect1[1] = z -TriVertex[2];
    r1[0]        = A[0][0]*tempVect1[0] +A[0][1]*tempVect1[1];
    r1[1]        = A[1][0]*tempVect1[0] +A[1][1]*tempVect1[1];
    y1           = r1[0];                   z1           = r1[1];
    // Transform the in-plane slip vector components from TDCS into ADCS
    tempVect1[0] = by;                      tempVect1[1] = bz;
    r2[0]        = A[0][0]*tempVect1[0] +A[0][1]*tempVect1[1];
    r2[1]        = A[1][0]*tempVect1[0] +A[1][1]*tempVect1[1];
    by1          = r2[0];                   bz1          = r2[1];
    
    // Calculate strains associated with an angular dislocation in ADCS
    AngDisStrain_inStrnHS(x, y1, z1, (-1.0*M_PI+alpha), bx, by1, bz1, nu, e); 
    // Transform strains from ADCS into TDCS
    B[0][0] = 1.0;          B[0][1] = 0.0;      B[0][2] = 0.0;
    B[1][0] = 0.0;          B[1][1] = A[0][0];  B[1][2] = A[1][0];
    B[2][0] = 0.0;          B[2][1] = A[0][1];  B[2][2] = A[1][1];// 3x3 Transformation matrix
    
    TensTrans_inStrnHS(e, e_out, B); //the e_out is then send back from the function

    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void AngDisStrain_inStrnHS(double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double e[6])
// AngDisStrain calculates the strains associated with an angular dislocation in an elastic full-space.
{   double       sinA,           cosA,           eta,            zeta;
    double       x2,             y2,             z2,             r2;
    double       r,              r3,             rz,             r2z2;
    double       r3z,            W,              W2,             Wr;
    double       W2r,            Wr3,            W2r2,           C;
    double       S,              rFi_rx,         rFi_ry,         rFi_rz;
    double       Exx,            Exy,            Exz,            Eyy;
    double       Eyz,            Ezz;
    //----------------------------------------------------------   

    sinA = sin(alpha);                          cosA = cos(alpha);
    eta = y*cosA - z*sinA;                      zeta = y*sinA + z*cosA;

    x2 = x*x;               y2   = y*y;                 z2  = z*z;
    r2 = x2 + y2 + z2;      r    = pow(r2,0.5);         r3  = r*r2;
    rz = r*(r-z);           r2z2 = r2*pow((r-z),2.0);   r3z = r3*(r-z);

    W  = zeta-r;            W2   = W*W;                 Wr  = W*r;
    W2r= W2*r;              Wr3  = W*r3;                W2r2= W2*r2;

    C = (r*cosA-z)/Wr;      S    = (r*sinA-y)/Wr;
    // Partial derivatives of the Burgers' function
    rFi_rx = (eta/r/(r-zeta)  -y/r/(r-z))   /4.0/M_PI;
    rFi_ry = (  x/r/(r-z)-cosA*x/r/(r-zeta))/4.0/M_PI;
    rFi_rz = (sinA*x/r/(r-zeta))/4.0/M_PI;
    //----------------------------------------------------------
    Exx = bx*(rFi_rx)  +bx/8.0/M_PI/(1.0-nu)*(eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-x2*y/r2z2-x2*y/r3z)
                       -by*x/8.0/M_PI/(1.0-nu)*(((2.0*nu+1.0)/Wr+x2/W2r2-x2/Wr3)*cosA+(2.0*nu+1.0)/rz-x2/r2z2-x2/r3z)
                       +bz*x*sinA/8.0/M_PI/(1.0-nu)*((2.0*nu+1.0)/Wr+x2/W2r2-x2/Wr3);

    Eyy = by*(rFi_ry)  +bx/8.0/M_PI/(1.0-nu)*((1.0/Wr+S*S-y2/Wr3)*eta+(2.0*nu+1.0)*y/rz-pow(y,3.0)/r2z2-pow(y,3.0)/r3z-2.0*nu*cosA*S)
                       -by*x/8.0/M_PI/(1.0-nu)*(1.0/rz-y2/r2z2-y2/r3z+(1.0/Wr+S*S-y2/Wr3)*cosA)
                       +bz*x*sinA/8.0/M_PI/(1.0-nu)*(1.0/Wr+S*S-y2/Wr3);

    Ezz = bz*(rFi_rz)+bx/8.0/M_PI/(1.0-nu)*(eta/W/r+eta*C*C-eta*z2/Wr3+y*z/r3+2.0*nu*sinA*C)
                      -by*x/8.0/M_PI/(1.0-nu)*((1.0/Wr+C*C-z2/Wr3)*cosA+z/r3)
                      +bz*x*sinA/8.0/M_PI/(1.0-nu)*(1.0/Wr+C*C-z2/Wr3);

    Exy = bx*(rFi_ry)/2.0+by*(rFi_rx)/2.0  -bx/8.0/M_PI/(1.0-nu)*(x*y2/r2z2-nu*x/rz+x*y2/r3z-nu*x*cosA/Wr+eta*x*S/Wr+eta*x*y/Wr3)
                                           +by/8.0/M_PI/(1.0-nu)*(x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+x2*y*cosA/Wr3+x2*cosA*S/Wr)
                                           -bz*sinA/8.0/M_PI/(1.0-nu)*(nu*S+x2*S/Wr+x2*y/Wr3);

    Exz = bx*(rFi_rz)/2.0+bz*(rFi_rx)/2.0  -bx/8.0/M_PI/(1.0-nu)*(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+eta*x*z/Wr3)
                                           +by/8.0/M_PI/(1.0-nu)*(-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+x2*cosA*C/Wr)
                                           -bz*sinA/8.0/M_PI/(1.0-nu)*(nu*C+x2*C/Wr+x2*z/Wr3);

    Eyz = by*(rFi_rz)/2.0+bz*(rFi_ry)/2.0  +bx/8.0/M_PI/(1.0-nu)*(y2/r3-nu/r-nu*cosA*C+nu*sinA*S+eta*sinA*cosA/W2-eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3)
                                           -by*x/8.0/M_PI/(1.0-nu)*(y/r3+sinA*cosA*cosA/W2-cosA*(y*cosA+z*sinA)/W2r+y*z*cosA/W2r2-y*z*cosA/Wr3)
                                           -bz*x*sinA/8.0/M_PI/(1.0-nu)*(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/W2r-y*z/W2r2);

    e[0] = Exx;             e[1] = Exy;             e[2] = Exz;
    e[3] = Eyy;             e[4] = Eyz;             e[5] = Ezz;

    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void TensTrans_inStrnHS(double e_in[6], double e_out[6], double B[3][3])
{
    // TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate system to x2y2z2 coordinate system. "A" is the transformation matrix, 
    // whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
    // of A (i.e., A') does the transformation from x2y2z2 into x1y1z1.
    double Txx1,         Txy1,           Txz1,           Tyy1;
    double Tyz1,         Tzz1,           Txx2,           Txy2;
    double Txz2,         Tyy2,           Tyz2,           Tzz2;
    double A[9];
    //---------------------------------
    Txx1 = e_in[0];         Txy1 = e_in[1];         Txz1 = e_in[2];
    Tyy1 = e_in[3];         Tyz1 = e_in[4];         Tzz1 = e_in[5];

    A[0] = B[0][0];         A[1] = B[1][0];         A[2] = B[2][0];
    A[3] = B[0][1];         A[4] = B[1][1];         A[5] = B[2][1];
    A[6] = B[0][2];         A[7] = B[1][2];         A[8] = B[2][2];

    Txx2 = A[0]*A[0]*Txx1 +         2.0*A[0]*A[3] *Txy1 +          2.0*A[0]*A[6] *Txz1 +          2.0*A[3]*A[6] *Tyz1 + A[3]*A[3]*Tyy1 + A[6]*A[6]*Tzz1;
    Tyy2 = A[1]*A[1]*Txx1 +         2.0*A[1]*A[4] *Txy1 +          2.0*A[1]*A[7] *Txz1 +          2.0*A[4]*A[7] *Tyz1 + A[4]*A[4]*Tyy1 + A[7]*A[7]*Tzz1;
    Tzz2 = A[2]*A[2]*Txx1 +         2.0*A[2]*A[5] *Txy1 +          2.0*A[2]*A[8] *Txz1 +          2.0*A[5]*A[8] *Tyz1 + A[5]*A[5]*Tyy1 + A[8]*A[8]*Tzz1;
    Txy2 = A[0]*A[1]*Txx1 + (A[0]*A[4]+ A[1]*A[3])*Txy1 + (A[0]*A[7] + A[1]*A[6])*Txz1 + (A[7]*A[3] + A[6]*A[4])*Tyz1 + A[4]*A[3]*Tyy1 + A[6]*A[7]*Tzz1;
    Txz2 = A[0]*A[2]*Txx1 + (A[0]*A[5]+ A[2]*A[3])*Txy1 + (A[0]*A[8] + A[2]*A[6])*Txz1 + (A[8]*A[3] + A[6]*A[5])*Tyz1 + A[5]*A[3]*Tyy1 + A[6]*A[8]*Tzz1;
    Tyz2 = A[1]*A[2]*Txx1 + (A[2]*A[4]+ A[1]*A[5])*Txy1 + (A[2]*A[7] + A[1]*A[8])*Txz1 + (A[7]*A[5] + A[8]*A[4])*Tyz1 + A[4]*A[5]*Tyy1 + A[7]*A[8]*Tzz1;
    //---------------------------------
    e_out[0] = Txx2;        e_out[1] = Txy2;      e_out[2] = Txz2;
    e_out[3] = Tyy2;        e_out[4] = Tyz2;      e_out[5] = Tzz2;
    
    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void    AngSetupFSC_S_inStrnHS(double Stress[6],double Strain[6], double X,double Y,double Z,double bX,double bY,double bZ,double PA[3], double PB[3],double mu,double lambda)
{   // AngSetupFSC_S calculates the Free Surface Correction to strains and  stresses associated with angular dislocation pair on each TD side.
    int   i,        I;
    double nu,      beta,   TempVal1,       TempVal2,       eps;
    
    double SideVec[3],       eZ[3],          ey1[3],         ey2[3],     ey3[3];
    double TempVect1[3],     TempVect2[3],   yA[3],          yB[3];
    double v1A[6],           v1B[6],         v_vals[6];

    double A[3][3],          A_t[3][3];
    eps = DBL_EPSILON;//maybe switch to doubles??
    nu  = lambda/(mu+lambda)/2.0; // Poisson's ratio
    //Calculate TD side vector and the angle of the angular dislocation pair
    SideVec[0] = PB[0]-PA[0];           SideVec[1] = PB[1]-PA[1];       SideVec[2] = PB[2]-PA[2];
    eZ[0]      = 0.0;                   eZ[1]      = 0.0;               eZ[2]      = 1.0;

    TempVal1   = pow( (pow(SideVec[0],2.0) + pow(SideVec[1],2.0)   + pow(SideVec[2],2.0)),0.5);
    TempVal2   =  -1.0*SideVec[0]*eZ[0]    + -1.0*SideVec[1]*eZ[1] + -1.0*SideVec[2]*eZ[2];    
    beta       = acos(TempVal2/TempVal1);
    if ((fabs(beta) < eps) || (fabs(M_PI-beta)< eps))
    {   for (i = 0; i < 6; i++)         {       Stress[i] = 0.0;        Strain[i] = 0.0;            } 
    }
    else
    {
        ey1[0]   = SideVec[0];          ey1[1]   = SideVec[1];           ey1[2]   = 0.0;
        TempVal1 = pow( (pow(ey1[0],2.0) + pow(ey1[1],2.0)   + pow(ey1[2],2.0)),0.5);
        ey1[0]  /= TempVal1;            ey1[1]  /= TempVal1;             ey1[2]  /= TempVal1;
        ey3[0]   = -1.0*eZ[0];          ey3[1]   = -1.0*eZ[1];          ey3[2]    = -1.0*eZ[2];
        
        ey2[0]   = ey3[1]*ey1[2] -ey3[2]*ey1[1]; 
        ey2[1]   = ey3[2]*ey1[0] -ey3[0]*ey1[2]; 
        ey2[2]   = ey3[0]*ey1[1] -ey3[1]*ey1[0]; 
        TempVal1 = pow( (pow(ey2[0],2.0) + pow(ey2[1],2.0)   + pow(ey2[2],2.0)),0.5);
        ey2[0]  /= TempVal1;            ey2[1]  /= TempVal1;             ey2[2]  /= TempVal1;
        //Transformation matrix
        A[0][0]   = ey1[0];              A[0][1]   = ey2[0];              A[0][2]   = ey3[0];
        A[1][0]   = ey1[1];              A[1][1]   = ey2[1];              A[1][2]   = ey3[1];
        A[2][0]   = ey1[2];              A[2][1]   = ey2[2];              A[2][2]   = ey3[2];
        
        A_t[0][0] = ey1[0];              A_t[0][1] = ey1[1];              A_t[0][2] = ey1[2];
        A_t[1][0] = ey2[0];              A_t[1][1] = ey2[1];              A_t[1][2] = ey2[2];
        A_t[2][0] = ey3[0];              A_t[2][1] = ey3[1];              A_t[2][2] = ey3[2];
       
        //Transform coordinates from EFCS to the first ADCS
        CoordTrans_inStrnHS(yA, (X-PA[0]), (Y-PA[1]), (Z-PA[2]), A); 
        //Transform coordinates from EFCS to the second ADCS
        CoordTrans_inStrnHS(TempVect1, SideVec[0], SideVec[1], SideVec[2], A); 
        yB[0] = yA[0] - TempVect1[0];
        yB[1] = yA[1] - TempVect1[1];
        yB[2] = yA[2] - TempVect1[2];
        //Transform slip vector components from EFCS to ADCS
        CoordTrans_inStrnHS(TempVect2, bX, bY, bZ, A); 
        //Determine the best arteact-free configuration for the calculation points near the free furface
        I = (beta*yA[0]) >=0.0 ? 1 : 0;
        //For singularities at surface
        for (i = 0; i < 6; i++)         {       v1A[i] = 0.0;       v1B[i] = 0.0;           }
        // Configuration I
        if (I == 1)
        {   
            AngDisStrainFSC_inStrnHS((-1.0*yA[0]),(-1.0*yA[1]), yA[2], (M_PI-beta),(-1.0*TempVect2[0]),(-1.0*TempVect2[1]), TempVect2[2], nu,(-1.0*PA[2]), v1A);
            v1A[2] *= -1.0; //this is the strain_xz component (strain_13)
            v1A[4] *= -1.0; //this is the strain_yz component (strain_23)
            AngDisStrainFSC_inStrnHS((-1.0*yB[0]),(-1.0*yB[1]), yB[2], (M_PI-beta),(-1.0*TempVect2[0]),(-1.0*TempVect2[1]), TempVect2[2], nu,(-1.0*PB[2]), v1B);
            v1B[2] *= -1.0; //this is the strain_xz component (strain_13)
            v1B[4] *= -1.0; //this is the strain_yz component (strain_23)
        }
        // Configuration II
        else if (I == 0)
        {   
            AngDisStrainFSC_inStrnHS(yA[0],yA[1], yA[2], beta,TempVect2[0],TempVect2[1], TempVect2[2], nu,(-1.0*PA[2]), v1A);
            AngDisStrainFSC_inStrnHS(yB[0],yB[1], yB[2], beta,TempVect2[0],TempVect2[1], TempVect2[2], nu,(-1.0*PB[2]), v1B);
        }
        //Calculate total Free Surface Correction to strains in ADCS
    	for (i = 0; i < 6; i++)         {           v_vals[i] = v1B[i] - v1A[i];            }
		
    	//Transform total Free Surface Correction to strains from ADCS to EFCS
    	TensTrans_inStrnHS(v_vals, Strain, A_t);
		
    	Stress[0] = 2.0*mu*Strain[0] +lambda*(Strain[0] +Strain[3] +Strain[5]); //sig_xx
   	 	Stress[3] = 2.0*mu*Strain[3] +lambda*(Strain[0] +Strain[3] +Strain[5]); //sig_yy
    	Stress[5] = 2.0*mu*Strain[5] +lambda*(Strain[0] +Strain[3] +Strain[5]); //sig_zz
    	Stress[1] = 2.0*mu*Strain[1]; //sig_xy
    	Stress[2] = 2.0*mu*Strain[2]; //sig_xz
    	Stress[4] = 2.0*mu*Strain[4]; //sig_yz 
    }
    return;
}

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

void  AngDisStrainFSC_inStrnHS(double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a, double Strain[6])
{// AngDisStrainFSC calculates the harmonic function contribution to the strains associated with an angular dislocation in an elastic half-space.
    
    double sinB,         cosB,           cotB,           y3b,        z1b;
    double z3b,          rb2,            rb,             W1,         W2;
    double W3,           W4,             W5,             W6,         W7;
    double W8,           W9,             N1,             rFib_ry2,   rFib_ry1;
    double rFib_ry3,     v11,            v12,            v13,        v22;
    double v23,          v33;
    
    
    sinB = sin(beta);                   cosB = cos(beta);               cotB = cosB/sinB;//cot(beta);
    y3b  = y3+2.0*a;                    z1b  = y1*cosB+y3b*sinB;        z3b  = -y1*sinB+y3b*cosB;
    rb2  = y1*y1 + y2*y2 + y3b*y3b;     rb   = pow(rb2,0.5);            N1   = 1.0-2.0*nu;

    W1   = rb*cosB+y3b;                 W2   = cosB   +a/rb;            W3   = cosB +y3b/rb;
    W4   = nu +a/rb;                    W5   = 2.0*nu +a/rb;            W6   = rb   +y3b;
    W7   = rb +z3b;                     W8   = y3 +a;                   W9   = 1.0  +a/rb/cosB;

    // Partial derivatives of the Burgers' function
    rFib_ry2 =      z1b/rb/(rb+z3b)      -y1/rb/(rb+y3b); // y2 = x in ADCS
    rFib_ry1 =       y2/rb/(rb+y3b) -cosB*y2/rb/(rb+z3b); // y1 =y in ADCS
    rFib_ry3 = -sinB*y2/rb/(rb+z3b);                      // y3 = z in ADCS
    
//-------------------------------//-------------------------------

v11 = b1*(0.25*((-2.0+2.0*nu)*N1*rFib_ry1*pow(cotB,2.0)-N1*y2/pow(W6,2.0)*((1.0-W5)*cotB-y1/W6*W4)/rb*y1+N1*y2/W6*(a/pow(rb,3.0)*y1*cotB-1.0/W6*W4+pow(y1,2.0)/pow(W6,2.0)*W4/rb+pow(y1,2.0)/W6*a/pow(rb,3.0))-N1*y2*cosB*cotB/pow(W7,2.0)*W2*(y1/rb-sinB)-N1*y2*cosB*cotB/W7*a/pow(rb,3.0)*y1-3*a*y2*W8*cotB/pow(rb,5.0)*y1-y2*W8/pow(rb,3.0)/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1-y2*W8/rb2/pow(W6,2.0)*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1+y2*W8/rb/W6*(1/W6*W5-pow(y1,2.0)/pow(W6,2.0)*W5/rb-pow(y1,2.0)/W6*a/pow(rb,3.0)+a/rb2-2.0*a*pow(y1,2.0)/pow(rb2,2.0))-y2*W8/pow(rb,3.0)/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*y1-y2*W8/rb/pow(W7,2.0)*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*(y1/rb-sinB)+y2*W8/rb/W7*(-cosB/pow(W7,2.0)*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)*(y1/rb-sinB)+cosB/W7*(1.0/rb*cosB*y1*(N1*cosB-a/rb)*cotB+W1*a/pow(rb,3.0)*y1*cotB+(2.0-2.0*nu)*(1.0/rb*sinB*y1-1.0)*cosB)+2*a*y3b*cosB*cotB/pow(rb2,2.0)*y1))/M_PI/(1.0-nu))+
      b2*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)+nu)/rb*y1/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0)*cosB*(y1/rb-sinB)/W7)-N1/pow(W6,2.0)*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+pow(y1,2.0)/W6*W4)/rb*y1+N1/W6*(-N1*cotB+a*cotB/rb-a*pow(y1,2.0)*cotB/pow(rb,3.0)+2*y1/W6*W4-pow(y1,3.0)/pow(W6,2.0)*W4/rb-pow(y1,3.0)/W6*a/pow(rb,3.0))+N1*cotB/pow(W7,2.0)*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*(y1/rb-sinB)-N1*cotB/W7*(pow(cosB,2.0)-a*(1.0/rb*sinB*y1-1.0)/rb/cosB+a*(rb*sinB-y1)/pow(rb,3.0)/cosB*y1)-a*W8*cotB/pow(rb,3.0)+3*a*pow(y1,2.0)*W8*cotB/pow(rb,5.0)-W8/pow(W6,2.0)*(2*nu+1.0/rb*(N1*y1*cotB+a)-pow(y1,2.0)/rb/W6*W5-a*pow(y1,2.0)/pow(rb,3.0))/rb*y1+W8/W6*(-1.0/pow(rb,3.0)*(N1*y1*cotB+a)*y1+1.0/rb*N1*cotB-2.0*y1/rb/W6*W5+pow(y1,3.0)/pow(rb,3.0)/W6*W5+pow(y1,3.0)/rb2/pow(W6,2.0)*W5+pow(y1,3.0)/pow(rb2,2.0)/W6*a-2.0*a/pow(rb,3.0)*y1+3*a*pow(y1,3.0)/pow(rb,5.0))-W8*cotB/pow(W7,2.0)*(-cosB*sinB+a*y1*y3b/pow(rb,3.0)/cosB+(rb*sinB-y1)/rb*((2.0-2.0*nu)*cosB-W1/W7*W9))*(y1/rb-sinB)+W8*cotB/W7*(a*y3b/pow(rb,3.0)/cosB-3*a*pow(y1,2.0)*y3b/pow(rb,5.0)/cosB+(1.0/rb*sinB*y1-1.0)/rb*((2.0-2.0*nu)*cosB-W1/W7*W9)-(rb*sinB-y1)/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7*W9)*y1+(rb*sinB-y1)/rb*(-1.0/rb*cosB*y1/W7*W9+W1/pow(W7,2.0)*W9*(y1/rb-sinB)+W1/W7*a/pow(rb,3.0)/cosB*y1)))/M_PI/(1.0-nu))+
      b3*(0.25*(N1*(-y2/pow(W6,2.0)*(1.0+a/rb)/rb*y1-y2/W6*a/pow(rb,3.0)*y1+y2*cosB/pow(W7,2.0)*W2*(y1/rb-sinB)+y2*cosB/W7*a/pow(rb,3.0)*y1)+y2*W8/pow(rb,3.0)*(a/rb2+1.0/W6)*y1-y2*W8/rb*(-2.0*a/pow(rb2,2.0)*y1-1.0/pow(W6,2.0)/rb*y1)-y2*W8*cosB/pow(rb,3.0)/W7*(W1/W7*W2+a*y3b/rb2)*y1-y2*W8*cosB/rb/pow(W7,2.0)*(W1/W7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(1.0/rb*cosB*y1/W7*W2-W1/pow(W7,2.0)*W2*(y1/rb-sinB)-W1/W7*a/pow(rb,3.0)*y1-2.0*a*y3b/pow(rb2,2.0)*y1))/M_PI/(1.0-nu));

v22 = b1*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)-nu)/rb*y2/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0-2.0*nu)*cosB/rb*y2/W7)+N1/pow(W6,2.0)*(y1*cotB*(1.0-W5)+nu*y3b-a+pow(y2,2.0)/W6*W4)/rb*y2-N1/W6*(a*y1*cotB/pow(rb,3.0)*y2+2*y2/W6*W4-pow(y2,3.0)/pow(W6,2.0)*W4/rb-pow(y2,3.0)/W6*a/pow(rb,3.0))+N1*z1b*cotB/pow(W7,2.0)*W2/rb*y2+N1*z1b*cotB/W7*a/pow(rb,3.0)*y2+3*a*y2*W8*cotB/pow(rb,5.0)*y1-W8/pow(W6,2.0)*(-2.0*nu+1.0/rb*(N1*y1*cotB-a)+pow(y2,2.0)/rb/W6*W5+a*pow(y2,2.0)/pow(rb,3.0))/rb*y2+W8/W6*(-1.0/pow(rb,3.0)*(N1*y1*cotB-a)*y2+2*y2/rb/W6*W5-pow(y2,3.0)/pow(rb,3.0)/W6*W5-pow(y2,3.0)/rb2/pow(W6,2.0)*W5-pow(y2,3.0)/pow(rb2,2.0)/W6*a+2*a/pow(rb,3.0)*y2-3*a*pow(y2,3.0)/pow(rb,5.0))-W8/pow(W7,2.0)*(pow(cosB,2.0)-1.0/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/pow(rb,3.0)-1.0/rb/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1))/rb*y2+W8/W7*(1/pow(rb,3.0)*(N1*z1b*cotB+a*cosB)*y2-3*a*y3b*z1b*cotB/pow(rb,5.0)*y2+1.0/pow(rb,3.0)/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*y2+1.0/rb2/pow(W7,2.0)*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*y2-1.0/rb/W7*(2*y2*pow(cosB,2.0)+a*z1b*cotB/pow(rb,3.0)*W1*y2-a*z1b*cotB/rb2*cosB*y2)))/M_PI/(1.0-nu))+
      b2*(0.25*((2.0-2.0*nu)*N1*rFib_ry2*pow(cotB,2.0)+N1/W6*((W5-1.0)*cotB+y1/W6*W4)-N1*pow(y2,2.0)/pow(W6,2.0)*((W5-1.0)*cotB+y1/W6*W4)/rb+N1*y2/W6*(-a/pow(rb,3.0)*y2*cotB-y1/pow(W6,2.0)*W4/rb*y2-y2/W6*a/pow(rb,3.0)*y1)-N1*cotB/W7*W9+N1*pow(y2,2.0)*cotB/pow(W7,2.0)*W9/rb+N1*pow(y2,2.0)*cotB/W7*a/pow(rb,3.0)/cosB-a*W8*cotB/pow(rb,3.0)+3*a*pow(y2,2.0)*W8*cotB/pow(rb,5.0)+W8/rb/W6*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))-pow(y2,2.0)*W8/pow(rb,3.0)/W6*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))-pow(y2,2.0)*W8/rb2/pow(W6,2.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))+y2*W8/rb/W6*(2*nu*y1/pow(W6,2.0)/rb*y2+a*y1/pow(rb,3.0)*(1.0/rb+1.0/W6)*y2-a*y1/rb*(-1.0/pow(rb,3.0)*y2-1.0/pow(W6,2.0)/rb*y2))+W8*cotB/rb/W7*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)-pow(y2,2.0)*W8*cotB/pow(rb,3.0)/W7*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)-pow(y2,2.0)*W8*cotB/rb2/pow(W7,2.0)*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)+y2*W8*cotB/rb/W7*(1.0/rb*cosB*y2/W7*W9-W1/pow(W7,2.0)*W9/rb*y2-W1/W7*a/pow(rb,3.0)/cosB*y2-2.0*a*y3b/pow(rb2,2.0)/cosB*y2))/M_PI/(1.0-nu))+
      b3*(0.25*(N1*(-sinB/rb*y2/W7+y2/pow(W6,2.0)*(1.0+a/rb)/rb*y1+y2/W6*a/pow(rb,3.0)*y1-z1b/pow(W7,2.0)*W2/rb*y2-z1b/W7*a/pow(rb,3.0)*y2)-y2*W8/pow(rb,3.0)*(a/rb2+1.0/W6)*y1+y1*W8/rb*(-2.0*a/pow(rb2,2.0)*y2-1.0/pow(W6,2.0)/rb*y2)+W8/pow(W7,2.0)*(sinB*(cosB-a/rb)+z1b/rb*(1.0+a*y3b/rb2)-1.0/rb/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1))/rb*y2-W8/W7*(sinB*a/pow(rb,3.0)*y2-z1b/pow(rb,3.0)*(1.0+a*y3b/rb2)*y2-2.0*z1b/pow(rb,5.0)*a*y3b*y2+1.0/pow(rb,3.0)/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*y2+1.0/rb2/pow(W7,2.0)*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*y2-1.0/rb/W7*(2*y2*cosB*sinB+a*z1b/pow(rb,3.0)*W1*y2-a*z1b/rb2*cosB*y2)))/M_PI/(1.0-nu));

v33 = b1*(0.25*((2.0-2.0*nu)*(N1*rFib_ry3*cotB-y2/pow(W6,2.0)*W5*(y3b/rb+1.0)-0.5*y2/W6*a/pow(rb,3.0)*2.0*y3b+y2*cosB/pow(W7,2.0)*W2*W3+0.5*y2*cosB/W7*a/pow(rb,3.0)*2.0*y3b)+y2/rb*(2*nu/W6+a/rb2)-0.5*y2*W8/pow(rb,3.0)*(2*nu/W6+a/rb2)*2.0*y3b+y2*W8/rb*(-2.0*nu/pow(W6,2.0)*(y3b/rb+1.0)-a/pow(rb2,2.0)*2.0*y3b)+y2*cosB/rb/W7*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)-0.5*y2*W8*cosB/pow(rb,3.0)/W7*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)*2.0*y3b-y2*W8*cosB/rb/pow(W7,2.0)*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)*W3+y2*W8*cosB/rb/W7*(-(cosB*y3b/rb+1.0)/W7*W2+W1/pow(W7,2.0)*W2*W3+0.5*W1/W7*a/pow(rb,3.0)*2.0*y3b-a/rb2+a*y3b/pow(rb2,2.0)*2.0*y3b))/M_PI/(1.0-nu))+
      b2*(0.25*((-2.0+2.0*nu)*N1*cotB*((y3b/rb+1.0)/W6-cosB*W3/W7)+(2.0-2.0*nu)*y1/pow(W6,2.0)*W5*(y3b/rb+1.0)+0.5*(2.0-2.0*nu)*y1/W6*a/pow(rb,3.0)*2.0*y3b+(2.0-2.0*nu)*sinB/W7*W2-(2.0-2.0*nu)*z1b/pow(W7,2.0)*W2*W3-0.5*(2.0-2.0*nu)*z1b/W7*a/pow(rb,3.0)*2.0*y3b+1.0/rb*(N1*cotB-2.0*nu*y1/W6-a*y1/rb2)-0.5*W8/pow(rb,3.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb2)*2.0*y3b+W8/rb*(2*nu*y1/pow(W6,2.0)*(y3b/rb+1.0)+a*y1/pow(rb2,2.0)*2.0*y3b)-1.0/W7*(cosB*sinB+W1*cotB/rb*((2.0-2.0*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))+W8/pow(W7,2.0)*(cosB*sinB+W1*cotB/rb*((2.0-2.0*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*W3-W8/W7*((cosB*y3b/rb+1.0)*cotB/rb*((2.0-2.0*nu)*cosB-W1/W7)-0.5*W1*cotB/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7)*2.0*y3b+W1*cotB/rb*(-(cosB*y3b/rb+1.0)/W7+W1/pow(W7,2.0)*W3)-0.5*a/pow(rb,3.0)*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*2.0*y3b+a/rb*(-z1b/rb2-y3b*sinB/rb2+y3b*z1b/pow(rb2,2.0)*2.0*y3b-sinB*W1/rb/W7-z1b*(cosB*y3b/rb+1.0)/rb/W7+0.5*z1b*W1/pow(rb,3.0)/W7*2.0*y3b+z1b*W1/rb/pow(W7,2.0)*W3)))/M_PI/(1.0-nu))+
      b3*(0.25*((2.0-2.0*nu)*rFib_ry3-(2.0-2.0*nu)*y2*sinB/pow(W7,2.0)*W2*W3-0.5*(2.0-2.0*nu)*y2*sinB/W7*a/pow(rb,3.0)*2.0*y3b+y2*sinB/rb/W7*(1.0+W1/W7*W2+a*y3b/rb2)-0.5*y2*W8*sinB/pow(rb,3.0)/W7*(1.0+W1/W7*W2+a*y3b/rb2)*2.0*y3b-y2*W8*sinB/rb/pow(W7,2.0)*(1.0+W1/W7*W2+a*y3b/rb2)*W3+y2*W8*sinB/rb/W7*((cosB*y3b/rb+1.0)/W7*W2-W1/pow(W7,2.0)*W2*W3-0.5*W1/W7*a/pow(rb,3.0)*2.0*y3b+a/rb2-a*y3b/pow(rb2,2.0)*2.0*y3b))/M_PI/(1.0-nu));

v12 = b1/2.0*(0.25*((-2.0+2.0*nu)*N1*rFib_ry2*pow(cotB,2.0)+N1/W6*((1.0-W5)*cotB-y1/W6*W4)-N1*pow(y2,2.0)/pow(W6,2.0)*((1.0-W5)*cotB-y1/W6*W4)/rb+N1*y2/W6*(a/pow(rb,3.0)*y2*cotB+y1/pow(W6,2.0)*W4/rb*y2+y2/W6*a/pow(rb,3.0)*y1)+N1*cosB*cotB/W7*W2-N1*pow(y2,2.0)*cosB*cotB/pow(W7,2.0)*W2/rb-N1*pow(y2,2.0)*cosB*cotB/W7*a/pow(rb,3.0)+a*W8*cotB/pow(rb,3.0)-3*a*pow(y2,2.0)*W8*cotB/pow(rb,5.0)+W8/rb/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-pow(y2,2.0)*W8/pow(rb,3.0)/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-pow(y2,2.0)*W8/rb2/pow(W6,2.0)*(-N1*cotB+y1/W6*W5+a*y1/rb2)+y2*W8/rb/W6*(-y1/pow(W6,2.0)*W5/rb*y2-y2/W6*a/pow(rb,3.0)*y1-2.0*a*y1/pow(rb2,2.0)*y2)+W8/rb/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-pow(y2,2.0)*W8/pow(rb,3.0)/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-pow(y2,2.0)*W8/rb2/pow(W7,2.0)*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)+y2*W8/rb/W7*(-cosB/pow(W7,2.0)*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)/rb*y2+cosB/W7*(1.0/rb*cosB*y2*(N1*cosB-a/rb)*cotB+W1*a/pow(rb,3.0)*y2*cotB+(2.0-2.0*nu)/rb*sinB*y2*cosB)+2*a*y3b*cosB*cotB/pow(rb2,2.0)*y2))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)+nu)/rb*y2/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0)*cosB/rb*y2/W7)-N1/pow(W6,2.0)*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+pow(y1,2.0)/W6*W4)/rb*y2+N1/W6*(-a*y1*cotB/pow(rb,3.0)*y2-pow(y1,2.0)/pow(W6,2.0)*W4/rb*y2-pow(y1,2.0)/W6*a/pow(rb,3.0)*y2)+N1*cotB/pow(W7,2.0)*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)/rb*y2-N1*cotB/W7*(-a/rb2*sinB*y2/cosB+a*(rb*sinB-y1)/pow(rb,3.0)/cosB*y2)+3*a*y2*W8*cotB/pow(rb,5.0)*y1-W8/pow(W6,2.0)*(2*nu+1.0/rb*(N1*y1*cotB+a)-pow(y1,2.0)/rb/W6*W5-a*pow(y1,2.0)/pow(rb,3.0))/rb*y2+W8/W6*(-1.0/pow(rb,3.0)*(N1*y1*cotB+a)*y2+pow(y1,2.0)/pow(rb,3.0)/W6*W5*y2+pow(y1,2.0)/rb2/pow(W6,2.0)*W5*y2+pow(y1,2.0)/pow(rb2,2.0)/W6*a*y2+3*a*pow(y1,2.0)/pow(rb,5.0)*y2)-W8*cotB/pow(W7,2.0)*(-cosB*sinB+a*y1*y3b/pow(rb,3.0)/cosB+(rb*sinB-y1)/rb*((2.0-2.0*nu)*cosB-W1/W7*W9))/rb*y2+W8*cotB/W7*(-3*a*y1*y3b/pow(rb,5.0)/cosB*y2+1.0/rb2*sinB*y2*((2.0-2.0*nu)*cosB-W1/W7*W9)-(rb*sinB-y1)/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7*W9)*y2+(rb*sinB-y1)/rb*(-1.0/rb*cosB*y2/W7*W9+W1/pow(W7,2.0)*W9/rb*y2+W1/W7*a/pow(rb,3.0)/cosB*y2)))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*(N1*(1/W6*(1.0+a/rb)-pow(y2,2.0)/pow(W6,2.0)*(1.0+a/rb)/rb-pow(y2,2.0)/W6*a/pow(rb,3.0)-cosB/W7*W2+pow(y2,2.0)*cosB/pow(W7,2.0)*W2/rb+pow(y2,2.0)*cosB/W7*a/pow(rb,3.0))-W8/rb*(a/rb2+1.0/W6)+pow(y2,2.0)*W8/pow(rb,3.0)*(a/rb2+1.0/W6)-y2*W8/rb*(-2.0*a/pow(rb2,2.0)*y2-1.0/pow(W6,2.0)/rb*y2)+W8*cosB/rb/W7*(W1/W7*W2+a*y3b/rb2)-pow(y2,2.0)*W8*cosB/pow(rb,3.0)/W7*(W1/W7*W2+a*y3b/rb2)-pow(y2,2.0)*W8*cosB/rb2/pow(W7,2.0)*(W1/W7*W2+a*y3b/rb2)+y2*W8*cosB/rb/W7*(1.0/rb*cosB*y2/W7*W2-W1/pow(W7,2.0)*W2/rb*y2-W1/W7*a/pow(rb,3.0)*y2-2.0*a*y3b/pow(rb2,2.0)*y2))/M_PI/(1.0-nu))+
      b1/2.0*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)-nu)/rb*y1/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0-2.0*nu)*cosB*(y1/rb-sinB)/W7)+N1/pow(W6,2.0)*(y1*cotB*(1.0-W5)+nu*y3b-a+pow(y2,2.0)/W6*W4)/rb*y1-N1/W6*((1.0-W5)*cotB+a*pow(y1,2.0)*cotB/pow(rb,3.0)-pow(y2,2.0)/pow(W6,2.0)*W4/rb*y1-pow(y2,2.0)/W6*a/pow(rb,3.0)*y1)-N1*cosB*cotB/W7*W2+N1*z1b*cotB/pow(W7,2.0)*W2*(y1/rb-sinB)+N1*z1b*cotB/W7*a/pow(rb,3.0)*y1-a*W8*cotB/pow(rb,3.0)+3*a*pow(y1,2.0)*W8*cotB/pow(rb,5.0)-W8/pow(W6,2.0)*(-2.0*nu+1.0/rb*(N1*y1*cotB-a)+pow(y2,2.0)/rb/W6*W5+a*pow(y2,2.0)/pow(rb,3.0))/rb*y1+W8/W6*(-1.0/pow(rb,3.0)*(N1*y1*cotB-a)*y1+1.0/rb*N1*cotB-pow(y2,2.0)/pow(rb,3.0)/W6*W5*y1-pow(y2,2.0)/rb2/pow(W6,2.0)*W5*y1-pow(y2,2.0)/pow(rb2,2.0)/W6*a*y1-3*a*pow(y2,2.0)/pow(rb,5.0)*y1)-W8/pow(W7,2.0)*(pow(cosB,2.0)-1.0/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/pow(rb,3.0)-1.0/rb/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1))*(y1/rb-sinB)+W8/W7*(1/pow(rb,3.0)*(N1*z1b*cotB+a*cosB)*y1-1.0/rb*N1*cosB*cotB+a*y3b*cosB*cotB/pow(rb,3.0)-3*a*y3b*z1b*cotB/pow(rb,5.0)*y1+1.0/pow(rb,3.0)/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*y1+1.0/rb/pow(W7,2.0)*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*(y1/rb-sinB)-1.0/rb/W7*(-a*cosB*cotB/rb*W1+a*z1b*cotB/pow(rb,3.0)*W1*y1-a*z1b*cotB/rb2*cosB*y1)))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*((2.0-2.0*nu)*N1*rFib_ry1*pow(cotB,2.0)-N1*y2/pow(W6,2.0)*((W5-1.0)*cotB+y1/W6*W4)/rb*y1+N1*y2/W6*(-a/pow(rb,3.0)*y1*cotB+1.0/W6*W4-pow(y1,2.0)/pow(W6,2.0)*W4/rb-pow(y1,2.0)/W6*a/pow(rb,3.0))+N1*y2*cotB/pow(W7,2.0)*W9*(y1/rb-sinB)+N1*y2*cotB/W7*a/pow(rb,3.0)/cosB*y1+3*a*y2*W8*cotB/pow(rb,5.0)*y1-y2*W8/pow(rb,3.0)/W6*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))*y1-y2*W8/rb2/pow(W6,2.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))*y1+y2*W8/rb/W6*(-2.0*nu/W6+2*nu*pow(y1,2.0)/pow(W6,2.0)/rb-a/rb*(1.0/rb+1.0/W6)+a*pow(y1,2.0)/pow(rb,3.0)*(1.0/rb+1.0/W6)-a*y1/rb*(-1.0/pow(rb,3.0)*y1-1.0/pow(W6,2.0)/rb*y1))-y2*W8*cotB/pow(rb,3.0)/W7*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*y1-y2*W8*cotB/rb/pow(W7,2.0)*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*(y1/rb-sinB)+y2*W8*cotB/rb/W7*(1.0/rb*cosB*y1/W7*W9-W1/pow(W7,2.0)*W9*(y1/rb-sinB)-W1/W7*a/pow(rb,3.0)/cosB*y1-2.0*a*y3b/pow(rb2,2.0)/cosB*y1))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*(N1*(-sinB*(y1/rb-sinB)/W7-1.0/W6*(1.0+a/rb)+pow(y1,2.0)/pow(W6,2.0)*(1.0+a/rb)/rb+pow(y1,2.0)/W6*a/pow(rb,3.0)+cosB/W7*W2-z1b/pow(W7,2.0)*W2*(y1/rb-sinB)-z1b/W7*a/pow(rb,3.0)*y1)+W8/rb*(a/rb2+1.0/W6)-pow(y1,2.0)*W8/pow(rb,3.0)*(a/rb2+1.0/W6)+y1*W8/rb*(-2.0*a/pow(rb2,2.0)*y1-1.0/pow(W6,2.0)/rb*y1)+W8/pow(W7,2.0)*(sinB*(cosB-a/rb)+z1b/rb*(1.0+a*y3b/rb2)-1.0/rb/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1))*(y1/rb-sinB)-W8/W7*(sinB*a/pow(rb,3.0)*y1+cosB/rb*(1.0+a*y3b/rb2)-z1b/pow(rb,3.0)*(1.0+a*y3b/rb2)*y1-2.0*z1b/pow(rb,5.0)*a*y3b*y1+1.0/pow(rb,3.0)/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*y1+1.0/rb/pow(W7,2.0)*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*(y1/rb-sinB)-1.0/rb/W7*(-a*cosB/rb*W1+a*z1b/pow(rb,3.0)*W1*y1-a*z1b/rb2*cosB*y1)))/M_PI/(1.0-nu));

v13 = b1/2.0*(0.25*((-2.0+2.0*nu)*N1*rFib_ry3*pow(cotB,2.0)-N1*y2/pow(W6,2.0)*((1.0-W5)*cotB-y1/W6*W4)*(y3b/rb+1.0)+N1*y2/W6*(0.5*a/pow(rb,3.0)*2.0*y3b*cotB+y1/pow(W6,2.0)*W4*(y3b/rb+1.0)+0.5*y1/W6*a/pow(rb,3.0)*2.0*y3b)-N1*y2*cosB*cotB/pow(W7,2.0)*W2*W3-0.5*N1*y2*cosB*cotB/W7*a/pow(rb,3.0)*2.0*y3b+a/pow(rb,3.0)*y2*cotB-1.5*a*y2*W8*cotB/pow(rb,5.0)*2.0*y3b+y2/rb/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-0.5*y2*W8/pow(rb,3.0)/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)*2.0*y3b-y2*W8/rb/pow(W6,2.0)*(-N1*cotB+y1/W6*W5+a*y1/rb2)*(y3b/rb+1.0)+y2*W8/rb/W6*(-y1/pow(W6,2.0)*W5*(y3b/rb+1.0)-0.5*y1/W6*a/pow(rb,3.0)*2.0*y3b-a*y1/pow(rb2,2.0)*2.0*y3b)+y2/rb/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-0.5*y2*W8/pow(rb,3.0)/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*2.0*y3b-y2*W8/rb/pow(W7,2.0)*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*W3+y2*W8/rb/W7*(-cosB/pow(W7,2.0)*(W1*(N1*cosB-a/rb)*cotB+(2.0-2.0*nu)*(rb*sinB-y1)*cosB)*W3+cosB/W7*((cosB*y3b/rb+1.0)*(N1*cosB-a/rb)*cotB+0.5*W1*a/pow(rb,3.0)*2.0*y3b*cotB+0.5*(2.0-2.0*nu)/rb*sinB*2.0*y3b*cosB)-a*cosB*cotB/rb2+a*y3b*cosB*cotB/pow(rb2,2.0)*2.0*y3b))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)+nu)*(y3b/rb+1.0)/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0)*cosB*W3/W7)-N1/pow(W6,2.0)*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+pow(y1,2.0)/W6*W4)*(y3b/rb+1.0)+N1/W6*(nu-0.5*a*y1*cotB/pow(rb,3.0)*2.0*y3b-pow(y1,2.0)/pow(W6,2.0)*W4*(y3b/rb+1.0)-0.5*pow(y1,2.0)/W6*a/pow(rb,3.0)*2.0*y3b)+N1*cotB/pow(W7,2.0)*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*W3-N1*cotB/W7*(cosB*sinB-0.5*a/rb2*sinB*2.0*y3b/cosB+0.5*a*(rb*sinB-y1)/pow(rb,3.0)/cosB*2.0*y3b)-a/pow(rb,3.0)*y1*cotB+1.5*a*y1*W8*cotB/pow(rb,5.0)*2.0*y3b+1.0/W6*(2.0*nu+1.0/rb*(N1*y1*cotB+a)-pow(y1,2.0)/rb/W6*W5-a*pow(y1,2.0)/pow(rb,3.0))-W8/pow(W6,2.0)*(2.0*nu+1.0/rb*(N1*y1*cotB+a)-pow(y1,2.0)/rb/W6*W5-a*pow(y1,2.0)/pow(rb,3.0))*(y3b/rb+1.0)+W8/W6*(-0.5/pow(rb,3.0)*(N1*y1*cotB+a)*2.0*y3b+0.5*pow(y1,2.0)/pow(rb,3.0)/W6*W5*2.0*y3b+pow(y1,2.0)/rb/pow(W6,2.0)*W5*(y3b/rb+1.0)+0.5*pow(y1,2.0)/pow(rb2,2.0)/W6*a*2.0*y3b+1.5*a*pow(y1,2.0)/pow(rb,5.0)*2.0*y3b)+cotB/W7*(-cosB*sinB+a*y1*y3b/pow(rb,3.0)/cosB+(rb*sinB-y1)/rb*((2.0-2.0*nu)*cosB-W1/W7*W9))-W8*cotB/pow(W7,2.0)*(-cosB*sinB+a*y1*y3b/pow(rb,3.0)/cosB+(rb*sinB-y1)/rb*((2.0-2.0*nu)*cosB-W1/W7*W9))*W3+W8*cotB/W7*(a/pow(rb,3.0)/cosB*y1-1.5*a*y1*y3b/pow(rb,5.0)/cosB*2.0*y3b+0.5/rb2*sinB*2.0*y3b*((2.0-2.0*nu)*cosB-W1/W7*W9)-0.5*(rb*sinB-y1)/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7*W9)*2.0*y3b+(rb*sinB-y1)/rb*(-(cosB*y3b/rb+1.0)/W7*W9+W1/pow(W7,2.0)*W9*W3+0.5*W1/W7*a/pow(rb,3.0)/cosB*2.0*y3b)))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*(N1*(-y2/pow(W6,2.0)*(1.0+a/rb)*(y3b/rb+1.0)-0.5*y2/W6*a/pow(rb,3.0)*2.0*y3b+y2*cosB/pow(W7,2.0)*W2*W3+0.5*y2*cosB/W7*a/pow(rb,3.0)*2.0*y3b)-y2/rb*(a/rb2+1.0/W6)+0.5*y2*W8/pow(rb,3.0)*(a/rb2+1.0/W6)*2.0*y3b-y2*W8/rb*(-a/pow(rb2,2.0)*2.0*y3b-1.0/pow(W6,2.0)*(y3b/rb+1.0))+y2*cosB/rb/W7*(W1/W7*W2+a*y3b/rb2)-0.5*y2*W8*cosB/pow(rb,3.0)/W7*(W1/W7*W2+a*y3b/rb2)*2.0*y3b-y2*W8*cosB/rb/pow(W7,2.0)*(W1/W7*W2+a*y3b/rb2)*W3+y2*W8*cosB/rb/W7*((cosB*y3b/rb+1.0)/W7*W2-W1/pow(W7,2.0)*W2*W3-0.5*W1/W7*a/pow(rb,3.0)*2.0*y3b+a/rb2-a*y3b/pow(rb2,2.0)*2.0*y3b))/M_PI/(1.0-nu))+
      b1/2.0*(0.25*((2.0-2.0*nu)*(N1*rFib_ry1*cotB-y1/pow(W6,2.0)*W5/rb*y2-y2/W6*a/pow(rb,3.0)*y1+y2*cosB/pow(W7,2.0)*W2*(y1/rb-sinB)+y2*cosB/W7*a/pow(rb,3.0)*y1)-y2*W8/pow(rb,3.0)*(2*nu/W6+a/rb2)*y1+y2*W8/rb*(-2.0*nu/pow(W6,2.0)/rb*y1-2.0*a/pow(rb2,2.0)*y1)-y2*W8*cosB/pow(rb,3.0)/W7*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)*y1-y2*W8*cosB/rb/pow(W7,2.0)*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(-1.0/rb*cosB*y1/W7*W2+W1/pow(W7,2.0)*W2*(y1/rb-sinB)+W1/W7*a/pow(rb,3.0)*y1+2*a*y3b/pow(rb2,2.0)*y1))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*((-2.0+2.0*nu)*N1*cotB*(1.0/rb*y1/W6-cosB*(y1/rb-sinB)/W7)-(2.0-2.0*nu)/W6*W5+(2.0-2.0*nu)*pow(y1,2.0)/pow(W6,2.0)*W5/rb+(2.0-2.0*nu)*pow(y1,2.0)/W6*a/pow(rb,3.0)+(2.0-2.0*nu)*cosB/W7*W2-(2.0-2.0*nu)*z1b/pow(W7,2.0)*W2*(y1/rb-sinB)-(2.0-2.0*nu)*z1b/W7*a/pow(rb,3.0)*y1-W8/pow(rb,3.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb2)*y1+W8/rb*(-2.0*nu/W6+2*nu*pow(y1,2.0)/pow(W6,2.0)/rb-a/rb2+2*a*pow(y1,2.0)/pow(rb2,2.0))+W8/pow(W7,2.0)*(cosB*sinB+W1*cotB/rb*((2.0-2.0*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*(y1/rb-sinB)-W8/W7*(1.0/rb2*cosB*y1*cotB*((2.0-2.0*nu)*cosB-W1/W7)-W1*cotB/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7)*y1+W1*cotB/rb*(-1.0/rb*cosB*y1/W7+W1/pow(W7,2.0)*(y1/rb-sinB))-a/pow(rb,3.0)*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*y1+a/rb*(-y3b*cosB/rb2+2*y3b*z1b/pow(rb2,2.0)*y1-cosB*W1/rb/W7-z1b/rb2*cosB*y1/W7+z1b*W1/pow(rb,3.0)/W7*y1+z1b*W1/rb/pow(W7,2.0)*(y1/rb-sinB))))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*((2.0-2.0*nu)*rFib_ry1-(2.0-2.0*nu)*y2*sinB/pow(W7,2.0)*W2*(y1/rb-sinB)-(2.0-2.0*nu)*y2*sinB/W7*a/pow(rb,3.0)*y1-y2*W8*sinB/pow(rb,3.0)/W7*(1.0+W1/W7*W2+a*y3b/rb2)*y1-y2*W8*sinB/rb/pow(W7,2.0)*(1.0+W1/W7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*sinB/rb/W7*(1.0/rb*cosB*y1/W7*W2-W1/pow(W7,2.0)*W2*(y1/rb-sinB)-W1/W7*a/pow(rb,3.0)*y1-2.0*a*y3b/pow(rb2,2.0)*y1))/M_PI/(1.0-nu));

v23 = b1/2.0*(0.25*(N1*(((2.0-2.0*nu)*pow(cotB,2.0)-nu)*(y3b/rb+1.0)/W6-((2.0-2.0*nu)*pow(cotB,2.0)+1.0-2.0*nu)*cosB*W3/W7)+N1/pow(W6,2.0)*(y1*cotB*(1.0-W5)+nu*y3b-a+pow(y2,2.0)/W6*W4)*(y3b/rb+1.0)-N1/W6*(0.5*a*y1*cotB/pow(rb,3.0)*2.0*y3b+nu-pow(y2,2.0)/pow(W6,2.0)*W4*(y3b/rb+1.0)-0.5*pow(y2,2.0)/W6*a/pow(rb,3.0)*2.0*y3b)-N1*sinB*cotB/W7*W2+N1*z1b*cotB/pow(W7,2.0)*W2*W3+0.5*N1*z1b*cotB/W7*a/pow(rb,3.0)*2.0*y3b-a/pow(rb,3.0)*y1*cotB+1.5*a*y1*W8*cotB/pow(rb,5.0)*2.0*y3b+1.0/W6*(-2.0*nu+1.0/rb*(N1*y1*cotB-a)+pow(y2,2.0)/rb/W6*W5+a*pow(y2,2.0)/pow(rb,3.0))-W8/pow(W6,2.0)*(-2.0*nu+1.0/rb*(N1*y1*cotB-a)+pow(y2,2.0)/rb/W6*W5+a*pow(y2,2.0)/pow(rb,3.0))*(y3b/rb+1.0)+W8/W6*(-0.5/pow(rb,3.0)*(N1*y1*cotB-a)*2.0*y3b-0.5*pow(y2,2.0)/pow(rb,3.0)/W6*W5*2.0*y3b-pow(y2,2.0)/rb/pow(W6,2.0)*W5*(y3b/rb+1.0)-0.5*pow(y2,2.0)/pow(rb2,2.0)/W6*a*2.0*y3b-1.5*a*pow(y2,2.0)/pow(rb,5.0)*2.0*y3b)+1.0/W7*(pow(cosB,2.0)-1.0/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/pow(rb,3.0)-1.0/rb/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1))-W8/pow(W7,2.0)*(pow(cosB,2.0)-1.0/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/pow(rb,3.0)-1.0/rb/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1))*W3+W8/W7*(0.5/pow(rb,3.0)*(N1*z1b*cotB+a*cosB)*2.0*y3b-1.0/rb*N1*sinB*cotB+a*z1b*cotB/pow(rb,3.0)+a*y3b*sinB*cotB/pow(rb,3.0)-1.5*a*y3b*z1b*cotB/pow(rb,5.0)*2.0*y3b+0.5/pow(rb,3.0)/W7*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*2.0*y3b+1.0/rb/pow(W7,2.0)*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*W1)*W3-1.0/rb/W7*(-a*sinB*cotB/rb*W1+0.5*a*z1b*cotB/pow(rb,3.0)*W1*2.0*y3b-a*z1b*cotB/rb*(cosB*y3b/rb+1.0))))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*((2.0-2.0*nu)*N1*rFib_ry3*pow(cotB,2.0)-N1*y2/pow(W6,2.0)*((W5-1.0)*cotB+y1/W6*W4)*(y3b/rb+1.0)+N1*y2/W6*(-0.5*a/pow(rb,3.0)*2.0*y3b*cotB-y1/pow(W6,2.0)*W4*(y3b/rb+1.0)-0.5*y1/W6*a/pow(rb,3.0)*2.0*y3b)+N1*y2*cotB/pow(W7,2.0)*W9*W3+0.5*N1*y2*cotB/W7*a/pow(rb,3.0)/cosB*2.0*y3b-a/pow(rb,3.0)*y2*cotB+1.5*a*y2*W8*cotB/pow(rb,5.0)*2.0*y3b+y2/rb/W6*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))-0.5*y2*W8/pow(rb,3.0)/W6*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))*2.0*y3b-y2*W8/rb/pow(W6,2.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb*(1.0/rb+1.0/W6))*(y3b/rb+1.0)+y2*W8/rb/W6*(2*nu*y1/pow(W6,2.0)*(y3b/rb+1.0)+0.5*a*y1/pow(rb,3.0)*(1.0/rb+1.0/W6)*2.0*y3b-a*y1/rb*(-0.5/pow(rb,3.0)*2.0*y3b-1.0/pow(W6,2.0)*(y3b/rb+1.0)))+y2*cotB/rb/W7*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)-0.5*y2*W8*cotB/pow(rb,3.0)/W7*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*2.0*y3b-y2*W8*cotB/rb/pow(W7,2.0)*((-2.0+2.0*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*W3+y2*W8*cotB/rb/W7*((cosB*y3b/rb+1.0)/W7*W9-W1/pow(W7,2.0)*W9*W3-0.5*W1/W7*a/pow(rb,3.0)/cosB*2.0*y3b+a/rb2/cosB-a*y3b/pow(rb2,2.0)/cosB*2.0*y3b))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*(N1*(-sinB*W3/W7+y1/pow(W6,2.0)*(1.0+a/rb)*(y3b/rb+1.0)+0.5*y1/W6*a/pow(rb,3.0)*2.0*y3b+sinB/W7*W2-z1b/pow(W7,2.0)*W2*W3-0.5*z1b/W7*a/pow(rb,3.0)*2.0*y3b)+y1/rb*(a/rb2+1.0/W6)-0.5*y1*W8/pow(rb,3.0)*(a/rb2+1.0/W6)*2.0*y3b+y1*W8/rb*(-a/pow(rb2,2.0)*2.0*y3b-1.0/pow(W6,2.0)*(y3b/rb+1.0))-1.0/W7*(sinB*(cosB-a/rb)+z1b/rb*(1.0+a*y3b/rb2)-1.0/rb/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1))+W8/pow(W7,2.0)*(sinB*(cosB-a/rb)+z1b/rb*(1.0+a*y3b/rb2)-1.0/rb/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1))*W3-W8/W7*(0.5*sinB*a/pow(rb,3.0)*2.0*y3b+sinB/rb*(1.0+a*y3b/rb2)-0.5*z1b/pow(rb,3.0)*(1.0+a*y3b/rb2)*2.0*y3b+z1b/rb*(a/rb2-a*y3b/pow(rb2,2.0)*2.0*y3b)+0.5/pow(rb,3.0)/W7*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*2.0*y3b+1.0/rb/pow(W7,2.0)*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*W1)*W3-1.0/rb/W7*(-a*sinB/rb*W1+0.5*a*z1b/pow(rb,3.0)*W1*2.0*y3b-a*z1b/rb*(cosB*y3b/rb+1.0))))/M_PI/(1.0-nu))+
      b1/2.0*(0.25*((2.0-2.0*nu)*(N1*rFib_ry2*cotB+1.0/W6*W5-pow(y2,2.0)/pow(W6,2.0)*W5/rb-pow(y2,2.0)/W6*a/pow(rb,3.0)-cosB/W7*W2+pow(y2,2.0)*cosB/pow(W7,2.0)*W2/rb+pow(y2,2.0)*cosB/W7*a/pow(rb,3.0))+W8/rb*(2*nu/W6+a/rb2)-pow(y2,2.0)*W8/pow(rb,3.0)*(2*nu/W6+a/rb2)+y2*W8/rb*(-2.0*nu/pow(W6,2.0)/rb*y2-2.0*a/pow(rb2,2.0)*y2)+W8*cosB/rb/W7*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)-pow(y2,2.0)*W8*cosB/pow(rb,3.0)/W7*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)-pow(y2,2.0)*W8*cosB/rb2/pow(W7,2.0)*(1.0-2.0*nu-W1/W7*W2-a*y3b/rb2)+y2*W8*cosB/rb/W7*(-1.0/rb*cosB*y2/W7*W2+W1/pow(W7,2.0)*W2/rb*y2+W1/W7*a/pow(rb,3.0)*y2+2*a*y3b/pow(rb2,2.0)*y2))/M_PI/(1.0-nu))+
      b2/2.0*(0.25*((-2.0+2.0*nu)*N1*cotB*(1.0/rb*y2/W6-cosB/rb*y2/W7)+(2.0-2.0*nu)*y1/pow(W6,2.0)*W5/rb*y2+(2.0-2.0*nu)*y1/W6*a/pow(rb,3.0)*y2-(2.0-2.0*nu)*z1b/pow(W7,2.0)*W2/rb*y2-(2.0-2.0*nu)*z1b/W7*a/pow(rb,3.0)*y2-W8/pow(rb,3.0)*(N1*cotB-2.0*nu*y1/W6-a*y1/rb2)*y2+W8/rb*(2*nu*y1/pow(W6,2.0)/rb*y2+2*a*y1/pow(rb2,2.0)*y2)+W8/pow(W7,2.0)*(cosB*sinB+W1*cotB/rb*((2.0-2.0*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))/rb*y2-W8/W7*(1.0/rb2*cosB*y2*cotB*((2.0-2.0*nu)*cosB-W1/W7)-W1*cotB/pow(rb,3.0)*((2.0-2.0*nu)*cosB-W1/W7)*y2+W1*cotB/rb*(-cosB/rb*y2/W7+W1/pow(W7,2.0)/rb*y2)-a/pow(rb,3.0)*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*y2+a/rb*(2.0*y3b*z1b/pow(rb2,2.0)*y2-z1b/rb2*cosB*y2/W7+z1b*W1/pow(rb,3.0)/W7*y2+z1b*W1/rb2/pow(W7,2.0)*y2)))/M_PI/(1.0-nu))+
      b3/2.0*(0.25*((2.0-2.0*nu)*rFib_ry2+(2.0-2.0*nu)*sinB/W7*W2-(2.0-2.0*nu)*pow(y2,2.0)*sinB/pow(W7,2.0)*W2/rb-(2.0-2.0*nu)*pow(y2,2.0)*sinB/W7*a/pow(rb,3.0)+W8*sinB/rb/W7*(1.0+W1/W7*W2+a*y3b/rb2)-pow(y2,2.0)*W8*sinB/pow(rb,3.0)/W7*(1.0+W1/W7*W2+a*y3b/rb2)-pow(y2,2.0)*W8*sinB/rb2/pow(W7,2.0)*(1.0+W1/W7*W2+a*y3b/rb2)+y2*W8*sinB/rb/W7*(1.0/rb*cosB*y2/W7*W2-W1/pow(W7,2.0)*W2/rb*y2-W1/W7*a/pow(rb,3.0)*y2-2.0*a*y3b/pow(rb2,2.0)*y2))/M_PI/(1.0-nu));
//-------------------------------//-------------------------------

    Strain[0] = v11;           Strain[1] = v12;         Strain[2] = v13;
    Strain[3] = v22;           Strain[4] = v23;         Strain[5] = v33;
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void    CrossProd_inStrainHalf(double v_in1[3], double v_in2[3], double v_out[3])
{    
    v_out[0]   = v_in1[1]*v_in2[2] - v_in1[2]*v_in2[1];
    v_out[1]   = v_in1[2]*v_in2[0] - v_in1[0]*v_in2[2];
    v_out[2]   = v_in1[0]*v_in2[1] - v_in1[1]*v_in2[0];

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void     NormalizVect_inStrainHalf(double v[3])
{
    double vNorm;
    vNorm  = pow( (pow(v[0],2.0) + pow(v[1],2.0)  + pow(v[2],2.0)),0.5);
    v[0]  /= vNorm;                   v[1]  /= vNorm;                    v[2]  /= vNorm;		

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
