// inialize the libraries
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <float.h>
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void        TDdispFS_inDispHS(double Displ[3],    double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double nu);   
void  TDdisp_HarFunc_inDispHS(double Displ[3],    double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double nu);   
void      CoordTrans_inDispHS(double newVal[3],  double x, double y, double z, double RotMat[3][3]); 
void   TriModeFind_inDispHS(int    TriMode[1],  double x, double y, double z, double p1a, double p1b, double p2a, double p2b, double p3a, double p3b);
void        TDSetupD_inDispHS(double DispVect[3], double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double TriVertex[3], double SideVec[3]);
void     AngSetupFSC_inDispHS(double DispVect[3], double x, double y, double z, double bx, double by, double bz, double PA[3], double PB[3], double nu);
void      AngDisDisp_inDispHS(double DispVect[3], double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu); 
void   AngDisDispFSC_inDispHS(double DispVect[3], double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a); 
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
 void dispinhalfspace_(double Displ[3], double X,double Y,double Z,double P1[3],double P2[3],double P3[3],double Ss,double Ds,double Ts,double nu)
// TDdispHS 
// Calculates displacements associated with a triangular dislocation in an 
// elastic half-space.
//
// TD: Triangular Dislocation
// EFCS: Earth-Fixed Coordinate System
// TDCS: Triangular Dislocation Coordinate System
// ADCS: Angular Dislocation Coordinate System
// 
// INPUTS
// X, Y and Z: 
// Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
// must have the same size. 
//
// P1,P2 and P3:
// Coordinates of TD vertices in EFCS. 
// 
// Ss, Ds and Ts:
// TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
//
// nu:
// Poisson's ratio.
//
// OUTPUTS
// ue, un and uv:
// Calculated displacement vector components in EFCS. ue, un and uv have
// the same unit as Ss, Ds and Ts in the inputs.
// 
// 
// Example: Calculate and plot the first component of displacement vector 
// on a regular grid.
// 
// [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,-5);
// [ue,un,uv] = TDdispHS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 -2],-1,2,3,.25);
// h = surf(X,Y,reshape(ue,size(X)),'edgecolor','none');
// view(2)
// axis equal
// axis tight
// set(gcf,'renderer','painters')

// Reference journal article: 
// Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
// artefact-free solution. 
// Submitted to Geophysical Journal International 

// Copyright (c) 2014 Mehdi Nikkhoo
// 
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, including 
// without limitation the rights to use, copy, modify, merge, publish, 
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the 
// following conditions:
// 
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

// I appreciate any comments or bug reports.

// Mehdi Nikkhoo
// created: 2013.1.24
// Last modified: 2014.7.30
// 
// VolcanoTectonics Research Group
// Section 2.1, Physics of Earthquakes and Volcanoes
// Department 2, Physics of the Earth
// Helmholtz Centre Potsdam
// German Research Centre for Geosciences (GFZ)
// 
// email: 
// mehdi.nikkhoo@gfz-potsdam.de 
// mehdi.nikkhoo@gmail.com
// translated into C by Olaf Zielke, Spring 2015
{

    double DispMS[3],           DispFSC[3],                 DispIS[3];
    double P1mirr[3],           P2mirr[3],                  P3mirr[3];
    //---------------------------------------------------------------
    //---------------------------------------------------------------   
    DispMS[0]  = 0.0;           DispMS[1]  = 0.0;           DispMS[2]  = 0.0;                
    DispFSC[0] = 0.0;           DispFSC[1] = 0.0;           DispFSC[2] = 0.0;     
    DispIS[0]  = 0.0;           DispIS[1]  = 0.0;           DispIS[2]  = 0.0;             

    P1mirr[0] = P1[0];          P2mirr[0] = P2[0];          P3mirr[0] = P3[0];
    P1mirr[1] = P1[1];          P2mirr[1] = P2[1];          P3mirr[1] = P3[1];
    P1mirr[2] = -1.0*P1[2];     P2mirr[2] = -1.0*P2[2];     P3mirr[2] = -1.0*P3[2];
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    // Calculate main dislocation contribution to displacements
    TDdispFS_inDispHS(       DispMS, X, Y,  Z, P1, P2, P3, Ss, Ds, Ts, nu);   
    // Calculate harmonic function contribution to displacements
    TDdisp_HarFunc_inDispHS(DispFSC, X, Y,  Z, P1, P2, P3, Ss, Ds, Ts, nu);  
    // Calculate image dislocation contribution to displacements
    TDdispFS_inDispHS(       DispIS, X,Y,Z,P1mirr,P2mirr,P3mirr,Ss,Ds,Ts,nu);
    
    if ((P1mirr[2] == 0.0) && (P2mirr[2] == 0.0) && (P3mirr[2] == 0.0)) 
    {   DispIS[2] = -1.0*DispIS[2];
    }
    // Calculate the complete displacement vector components in EFCS
    Displ[0]  = DispMS[0] + DispFSC[0] + DispIS[0];
    Displ[1]  = DispMS[1] + DispFSC[1] + DispIS[1];
    Displ[2]  = DispMS[2] + DispFSC[2] + DispIS[2];

    if ((P1mirr[2] == 0.0) && (P2mirr[2] == 0.0) && (P3mirr[2] == 0.0)) 
    {   Displ[0] = -1.0*Displ[0];
        Displ[1] = -1.0*Displ[1];
        Displ[2] = -1.0*Displ[2];
    }

    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void TDdispFS_inDispHS(double DisplVect[3],    double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double by, double bz, double bx, double nu)  
{   // TDdispFS:  Calculates displacements associated with a triangular dislocation in an elastic full-space.
    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward, respectively.
    int    casepLog,        casenLog;
    
    double Tempdouble1,     A,                  B,                  C,                      na;
    double x,               y,                  z,                  nb,                     nc;
    double Tempdouble2,     Fi;
    
    
    int    TriMode[1];
    double Vnorm[3],        Vstrike[3],         Vdip[3],            tempVect1[3],           tempVect2[3];
    double eY[3],           eZ[3],              e12[3],             e13[3],                 e23[3];
    double p1[3],           p2[3],              p3[3],              d_vals1[3],             d_vals2[3];
    double d_vals3[3],      d_comb[3],          a[3],               b[3],                   c[3];
    
    double Amat[3][3];
    
    p1[0] = 0.0;            p1[1] = 0.0;        p1[2] = 0.0;
    p2[0] = 0.0;            p2[1] = 0.0;        p2[2] = 0.0;
    p3[0] = 0.0;            p3[1] = 0.0;        p3[2] = 0.0;
    eY[0] = 0.0;            eY[1] = 1.0;        eY[2] = 0.0;
    eZ[0] = 0.0;            eZ[1] = 0.0;        eZ[2] = 1.0;
    
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    tempVect1[0] = P2[0] -P1[0];                tempVect1[1] = P2[1] -P1[1];                tempVect1[2] = P2[2] -P1[2];
    tempVect2[0] = P3[0] -P1[0];                tempVect2[1] = P3[1] -P1[1];                tempVect2[2] = P3[2] -P1[2];
        
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
    
    Tempdouble1  = pow (  ( pow(Vnorm[0],2.0) + pow(Vnorm[1],2.0) +pow(Vnorm[2],2.0) ),0.5);
    Vnorm[0]     = Vnorm[0]/Tempdouble1;        Vnorm[1]     = Vnorm[1]/Tempdouble1;        Vnorm[2]     = Vnorm[2]/Tempdouble1;
   
    Vstrike[0]   = eZ[1]*Vnorm[2] - eZ[2]*Vnorm[1];
    Vstrike[1]   = eZ[2]*Vnorm[0] - eZ[0]*Vnorm[2];
    Vstrike[2]   = eZ[0]*Vnorm[1] - eZ[1]*Vnorm[0];
    //-------------------------
    Tempdouble1  = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    
    if (Tempdouble1 == 0.0)
    {   Vstrike[0] = eY[0]*Vnorm[2];            Vstrike[1] = eY[1]*Vnorm[2];                Vstrike[2] = eY[2]*Vnorm[2];        
        // For horizontal elements in case of half-space calculation!!! => Correct the strike vector of image dislocation only
        if (P1[2] > 0.0)
        {   Vstrike[0] = -1.0*Vstrike[0];       Vstrike[1] = -1.0*Vstrike[1];               Vstrike[2] = -1.0*Vstrike[2];
    }   }
    
    Tempdouble1= pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    Vstrike[0] = Vstrike[0]/Tempdouble1;        Vstrike[1] = Vstrike[1]/Tempdouble1;        Vstrike[2] = Vstrike[2]/Tempdouble1;
    
    Vdip[0]    = Vnorm[1]*Vstrike[2] - Vnorm[2]*Vstrike[1];
    Vdip[1]    = Vnorm[2]*Vstrike[0] - Vnorm[0]*Vstrike[2];
    Vdip[2]    = Vnorm[0]*Vstrike[1] - Vnorm[1]*Vstrike[0];
    
    //Transform coordinates and slip vector components from EFCS into TDCS
    Amat[0][0]    = Vnorm[0];                   Amat[0][1]    = Vnorm[1];                   Amat[0][2]    = Vnorm[2];
    Amat[1][0]    = Vstrike[0];                 Amat[1][1]    = Vstrike[1];                 Amat[1][2]    = Vstrike[2];
    Amat[2][0]    = Vdip[0];                    Amat[2][1]    = Vdip[1];                    Amat[2][2]    = Vdip[2];
    //-------------------------------------------------------------------------
    CoordTrans_inDispHS(tempVect1,(X-P2[0]),(Y-P2[1]),(Z-P2[2]),Amat);
    x             = tempVect1[0];               y             = tempVect1[1];               z             = tempVect1[2];
    CoordTrans_inDispHS(tempVect1,(P1[0]-P2[0]),(P1[1]-P2[1]),(P1[2]-P2[2]),Amat);
    p1[0]         = tempVect1[0];               p1[1]         = tempVect1[1];               p1[2]         = tempVect1[2];
    CoordTrans_inDispHS(tempVect1,(P3[0]-P2[0]),(P3[1]-P2[1]),(P3[2]-P2[2]),Amat);
    p3[0]         = tempVect1[0];               p3[1]         = tempVect1[1];               p3[2]         = tempVect1[2];     
    // Calculate the unit vectors along TD sides in TDCS
    Tempdouble1 = pow (  ( pow((p2[0]-p1[0]),2.0) + pow((p2[1]-p1[1]),2.0) +pow((p2[2]-p1[2]),2.0 )),0.5);
    e12[0]      = (p2[0]-p1[0])/Tempdouble1;     e12[1]        = (p2[1]-p1[1])/Tempdouble1;   e12[2]        = (p2[2]-p1[2])/Tempdouble1;     
    Tempdouble1 = pow (  ( pow((p3[0]-p1[0]),2.0) + pow((p3[1]-p1[1]),2.0) +pow((p3[2]-p1[2]),2.0 )),0.5);
    e13[0]      = (p3[0]-p1[0])/Tempdouble1;     e13[1]        = (p3[1]-p1[1])/Tempdouble1;   e13[2]        = (p3[2]-p1[2])/Tempdouble1;
    Tempdouble1 = pow (  ( pow((p3[0]-p2[0]),2.0) + pow((p3[1]-p2[1]),2.0) +pow((p3[2]-p2[2]),2.0 )),0.5);
    e23[0]      = (p3[0]-p2[0])/Tempdouble1;     e23[1]        = (p3[1]-p2[1])/Tempdouble1;   e23[2]        = (p3[2]-p2[2])/Tempdouble1; 
    // Calculate the TD angles    
    Tempdouble1 = e12[0]*e13[0] +e12[1]*e13[1] +e12[2]*e13[2];
    A           = acos(Tempdouble1);
    Tempdouble1 = -1.0*e12[0]*e23[0] + -1.0*e12[1]*e23[1] + -1.0*e12[2]*e23[2];
    B           = acos(Tempdouble1);
    Tempdouble1 = e23[0]*e13[0] +e23[1]*e13[1] +e23[2]*e13[2];
    C           = acos(Tempdouble1);
    //-------------------------------------------------------------------------
    // Determine the best arteact-free configuration for each calculation point
    ///////////////////////////////////////////////////////////////////////////////////////////// 
    TriModeFind_inDispHS(TriMode,y,z,x, p1[1], p1[2], p2[1], p2[2], p3[1], p3[2]);
    ///////////////////////////////////////////////////////////////////////////////////////////// 
    if (TriMode[0] == 1)       {       casepLog = 1;   casenLog = 0;        } 
    if (TriMode[0] ==-1)       {       casepLog = 0;   casenLog = 1;        } 
    if (TriMode[0] == 0)       {       casepLog = 0;   casenLog = 0;        } 
    if (casepLog == 1) // Configuration I
    {   // Calculate first angular dislocation contribution
        tempVect1[0] = -1.0*e13[0];         tempVect1[1] = -1.0*e13[1];         tempVect1[2] = -1.0*e13[2];              
        TDSetupD_inDispHS(d_vals1, x, y, z, A, bx, by, bz, nu, p1, tempVect1);
        // Calculate second angular dislocation contribution
        TDSetupD_inDispHS(d_vals2, x, y, z, B, bx, by, bz, nu, p2, e12);
        // Calculate third angular dislocation contribution
        TDSetupD_inDispHS(d_vals3, x, y, z, C, bx, by, bz, nu, p3, e23);  
    }
    if (casenLog == 1) // Configuration II
    {  // Calculate first angular dislocation contribution             
        TDSetupD_inDispHS(d_vals1, x, y, z, A, bx, by, bz, nu, p1, e13);
        // Calculate second angular dislocation contribution
        tempVect1[0] = -1.0*e12[0];         tempVect1[1] = -1.0*e12[1];         tempVect1[2] = -1.0*e12[2];      
        TDSetupD_inDispHS(d_vals2, x, y, z, B, bx, by, bz, nu, p2, tempVect1);
        // Calculate third angular dislocation contribution
        tempVect1[0] = -1.0*e23[0];         tempVect1[1] = -1.0*e23[1];         tempVect1[2] = -1.0*e23[2]; 
        TDSetupD_inDispHS(d_vals3, x, y, z, C, bx, by, bz, nu, p3, tempVect1);     
    }
    //-------------------------------------------------------------------------
    // Calculate the "incomplete" displacement vector components in TDCS
    if ((casenLog == 1) || (casepLog == 1))
    {   d_comb[0]       = d_vals1[0]+d_vals2[0]+d_vals3[0]; // slip in x
        d_comb[1]       = d_vals1[1]+d_vals2[1]+d_vals3[1]; // slip in y
        d_comb[2]       = d_vals1[2]+d_vals2[2]+d_vals3[2]; // slip in z
    }
    else
    {   d_comb[0]       = NAN; // slip in x => supposed to be "NaN" => have to check
        d_comb[1]       = NAN; // slip in y
        d_comb[2]       = NAN; // slip in z
    }
    //-------------------------------------------------------------------------
    // Calculate the Burgers' function contribution corresponding to the TD
    a[0]        = -x;                   a[1]        = (p1[1]-y);            a[2]        = (p1[2]-z);
    b[0]        = -x;                   b[1]        = -y;                   b[2]        = -z;
    c[0]        = -x;                   c[1]        = (p3[1]-y);            c[2]        = (p3[2]-z);
      
    na          = pow( (pow(a[0],2.0) + pow(a[1],2.0) + pow(a[2],2.0)  ),0.5);
    nb          = pow( (pow(b[0],2.0) + pow(b[1],2.0) + pow(b[2],2.0)  ),0.5);
    nc          = pow( (pow(c[0],2.0) + pow(c[1],2.0) + pow(c[2],2.0)  ),0.5);
    //---------------------------------
    Tempdouble1 = (a[0]*(b[1]*c[2] -b[2]*c[1]) - a[1]*(b[0]*c[2] -b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]));
    Tempdouble2 = (na*nb*nc + (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])*nc  + (a[0]*c[0]+a[1]*c[1]+a[2]*c[2])*nb + (b[0]*c[0]+b[1]*c[1]+b[2]*c[2])*na);
    Fi          = -2.0*atan2(Tempdouble1,Tempdouble2)/4.0/M_PI;
    //---------------------------------
    // Calculate the complete displacement vector components in TDCS
    d_comb[0]  += bx*Fi;               d_comb[1]  += by*Fi;               d_comb[2]  += bz*Fi;
    // Transform the complete displacement vector components from TDCS into EFCS
    Amat[0][0] = Vnorm[0];          Amat[0][1] = Vstrike[0];        Amat[0][2] = Vdip[0];
    Amat[1][0] = Vnorm[1];          Amat[1][1] = Vstrike[1];        Amat[1][2] = Vdip[1];
    Amat[2][0] = Vnorm[2];          Amat[2][1] = Vstrike[2];        Amat[2][2] = Vdip[2];
    //---------------------------------
    CoordTrans_inDispHS(tempVect1,d_comb[0],d_comb[1],d_comb[2],Amat);
    //---------------------------------
    DisplVect[0] = tempVect1[0];    DisplVect[1] = tempVect1[1];    DisplVect[2] = tempVect1[2];
    //---------------------------------
    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void TDdisp_HarFunc_inDispHS(double DisplVect[3], double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double by, double bz, double bx, double nu)  
{   // TDdisp_HarFunc calculates the harmonic function contribution to the displacements associated with a triangular dislocation in a half-space.
    // The function cancels the surface normal tractions induced by the main and image dislocations.

    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD  as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward,  respectively.

    double Tempdouble1;

    double Vnorm[3],           Vstrike[3],      Vdip[3],        eZ[3];
    double tempVect1[3],       tempVect2[3],    eY[3];
    double d_vals1[3],         d_vals2[3],      d_vals3[3];

    double Amat[3][3];
    //--------------------------------------
    eY[0] = 0.0;            eY[1] = 1.0;        eY[2] = 0.0;
    eZ[0] = 0.0;            eZ[1] = 0.0;        eZ[2] = 1.0;
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    tempVect1[0] = P2[0] -P1[0];                tempVect1[1] = P2[1] -P1[1];                tempVect1[2] = P2[2] -P1[2];
    tempVect2[0] = P3[0] -P1[0];                tempVect2[1] = P3[1] -P1[1];                tempVect2[2] = P3[2] -P1[2];
        
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
    
    Tempdouble1  = pow (  ( pow(Vnorm[0],2.0) + pow(Vnorm[1],2.0) +pow(Vnorm[2],2.0) ),0.5);
    Vnorm[0]     = Vnorm[0]/Tempdouble1;        Vnorm[1]     = Vnorm[1]/Tempdouble1;        Vnorm[2]     = Vnorm[2]/Tempdouble1;
   
    Vstrike[0]   = eZ[1]*Vnorm[2] - eZ[2]*Vnorm[1];
    Vstrike[1]   = eZ[2]*Vnorm[0] - eZ[0]*Vnorm[2];
    Vstrike[2]   = eZ[0]*Vnorm[1] - eZ[1]*Vnorm[0];
    //-------------------------
    Tempdouble1  = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    
    if (Tempdouble1 == 0.0)
    {   Vstrike[0] = eY[0]*Vnorm[2];            Vstrike[1] = eY[1]*Vnorm[2];                Vstrike[2] = eY[2]*Vnorm[2];          
    }   
    Tempdouble1  = pow (  ( pow(Vstrike[0],2.0) + pow(Vstrike[1],2.0) +pow(Vstrike[2],2.0) ),0.5);
    
    Vstrike[0]   = Vstrike[0]/Tempdouble1;        Vstrike[1] = Vstrike[1]/Tempdouble1;        Vstrike[2] = Vstrike[2]/Tempdouble1;
    
    Vdip[0]      = Vnorm[1]*Vstrike[2] - Vnorm[2]*Vstrike[1];
    Vdip[1]      = Vnorm[2]*Vstrike[0] - Vnorm[0]*Vstrike[2];
    Vdip[2]      = Vnorm[0]*Vstrike[1] - Vnorm[1]*Vstrike[0];
    // Transform slip vector components from TDCS into EFCS
    Amat[0][0] = Vnorm[0];          Amat[0][1] = Vstrike[0];        Amat[0][2] = Vdip[0];
    Amat[1][0] = Vnorm[1];          Amat[1][1] = Vstrike[1];        Amat[1][2] = Vdip[1];
    Amat[2][0] = Vnorm[2];          Amat[2][1] = Vstrike[2];        Amat[2][2] = Vdip[2];
    
    CoordTrans_inDispHS(tempVect1,bx,by,bz,Amat);
    // Calculate contribution of angular dislocation pair on each TD side 
    AngSetupFSC_inDispHS(d_vals1,X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P1,P2,nu); // Side P1P2
    AngSetupFSC_inDispHS(d_vals2,X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P2,P3,nu); // Side P2P3
    AngSetupFSC_inDispHS(d_vals3,X,Y,Z,tempVect1[0],tempVect1[1],tempVect1[2],P3,P1,nu); // Side P3P1

    // Calculate total harmonic function contribution to displacements
    DisplVect[0]  = d_vals1[0] + d_vals2[0] + d_vals3[0];
    DisplVect[1]  = d_vals1[1] + d_vals2[1] + d_vals3[1];
    DisplVect[2]  = d_vals1[2] + d_vals2[2] + d_vals3[2];
    
    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void CoordTrans_inDispHS(double newVal[3], double x1, double x2, double x3, double A[3][3])
{
    // CoordTrans transforms the coordinates of the vectors, from x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
    // transformation matrix, whose columns e1,e2 and e3 are the unit base vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
    // in X1X2X3. The transpose of A (i.e., A') will transform the coordinates  from X1X2X3 into x1x2x3.
    
    newVal[0] = A[0][0]*x1 + A[0][1]*x2 + A[0][2]*x3;     newVal[1] = A[1][0]*x1 + A[1][1]*x2 + A[1][2]*x3;       newVal[2] = A[2][0]*x1 + A[2][1]*x2 + A[2][2]*x3;

    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void TriModeFind_inDispHS(int TriMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b)
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

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void TDSetupD_inDispHS(double DispVect[3], double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double TriVertex[3], double SideVec[3])
{   // TDESetupD transforms coordinates of the calculation points as well as  slip vector components from ADCS into TDCS. It then calculates the  displacements in ADCS and transforms them into TDCS.
    double A[2][2];
    double r1[2];
    double r2[2];
    double r3[2];
    double y1;
    double z1;
    double by1;
    double bz1;
    double tempVect1[2];
    double disp[3];
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
    // Calculate displacements associated with an angular dislocation in ADCS
    AngDisDisp_inDispHS(disp, x, y1, z1 , (-M_PI+alpha),bx, by1, bz1, nu); 
    // Transform displacements from ADCS into TDCS

    r3[0]        = A[0][0]*disp[1] +A[1][0]*disp[2]; // there the transposed of A is used!
    r3[1]        = A[0][1]*disp[1] +A[1][1]*disp[2];
    
    DispVect[0]  = disp[0];
    DispVect[1]  = r3[0];
    DispVect[2]  = r3[1];

    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void     AngSetupFSC_inDispHS(double DispVect[3], double x, double y, double z, double bx, double by, double bz, double PA[3], double PB[3], double nu)
{   // AngSetupFSC calculates the Free Surface Correction to displacements associated with angular dislocation pair on each TD side.
    int    i,                   I;
    double beta,                TempVal1,       TempVal2,           eps;

    double SideVec[3],          eZ[3],          d_vals[3],          TempVect1[3];
    double ey1[3],              ey2[3],         ey3[3],             yA[3];
    double TempVect2[3],        yB[3],          vA[3],              vB[3];
    
    double Amat[3][3],          Amat_t[3][3];
    
    
    eps = DBL_EPSILON;//maybe switch to doubles??
    //-------------------------------
    SideVec[0] = PB[0]-PA[0];           SideVec[1] = PB[1]-PA[1];       SideVec[2] = PB[2]-PA[2];
    eZ[0]      = 0.0;                   eZ[1]      = 0.0;               eZ[2]      = 1.0;

    // Calculate TD side vector and the angle of the angular dislocation pair
    TempVal1   = pow( (pow(SideVec[0],2.0) + pow(SideVec[1],2.0)   + pow(SideVec[2],2.0)),0.5);
    TempVal2   =  -1.0*SideVec[0]*eZ[0]    + -1.0*SideVec[1]*eZ[1] + -1.0*SideVec[2]*eZ[2];    
    beta       = acos(TempVal2/TempVal1);

    if ((fabs(beta) < eps) || (fabs(M_PI-beta)< eps))
    {   for (i = 0; i < 3; i++)         {      DispVect[i] = 0.0;            } 
    }
    else
    {
        ey1[0]   = SideVec[0];          ey1[1]   = SideVec[1];           ey1[2]   = 0.0;
        TempVal1 = pow( (pow(ey1[0],2.0) + pow(ey1[1],2.0)   + pow(ey1[2],2.0)),0.5);
        ey1[0]  /= TempVal1;            ey1[1]  /= TempVal1;             ey1[2]  /= TempVal1;
        ey3[0]   = -1.0*eZ[0];          ey3[1]   = -1.0*eZ[1];           ey3[2]   = -1.0*eZ[2];
        
        ey2[0]   = ey3[1]*ey1[2] -ey3[2]*ey1[1]; 
        ey2[1]   = ey3[2]*ey1[0] -ey3[0]*ey1[2]; 
        ey2[2]   = ey3[0]*ey1[1] -ey3[1]*ey1[0]; 
        //Transformation matrix
        Amat[0][0]   = ey1[0];              Amat[0][1]   = ey2[0];              Amat[0][2]   = ey3[0];
        Amat[1][0]   = ey1[1];              Amat[1][1]   = ey2[1];              Amat[1][2]   = ey3[1];
        Amat[2][0]   = ey1[2];              Amat[2][1]   = ey2[2];              Amat[2][2]   = ey3[2];
        
        Amat_t[0][0] = ey1[0];              Amat_t[0][1] = ey1[1];              Amat_t[0][2] = ey1[2];
        Amat_t[1][0] = ey2[0];              Amat_t[1][1] = ey2[1];              Amat_t[1][2] = ey2[2];
        Amat_t[2][0] = ey3[0];              Amat_t[2][1] = ey3[1];              Amat_t[2][2] = ey3[2];
        // Transform coordinates from EFCS to the first ADCS
        CoordTrans_inDispHS(yA, (x-PA[0]),(y-PA[1]),(z-PA[2]),Amat);
        // Transform coordinates from EFCS to the second ADCS
        CoordTrans_inDispHS(TempVect1, SideVec[0], SideVec[1], SideVec[2], Amat);
        yB[0]       = yA[0] -TempVect1[0];  yB[1]       = yA[1] -TempVect1[1];  yB[2]       = yA[2] -TempVect1[2]; 
    
        // Transform slip vector components from EFCS to ADCS
        CoordTrans_inDispHS(TempVect2, bx, by, bz, Amat);
        
        //Determine the best arteact-free configuration for the calculation points near the free furface
        I = (beta*yA[0]) >=0.0 ? 1 : 0;
        //For singularities at surface
        for (i = 0; i < 3; i++)         {       vA[i] = 0.0;       vB[i] = 0.0;           }

        if (I == 1)   // Configuration I
        {   AngDisDispFSC_inDispHS(vA, yA[0], yA[1], yA[2], (-M_PI+beta), TempVect2[0], TempVect2[1], TempVect2[2], nu, -PA[2]);
            AngDisDispFSC_inDispHS(vB, yB[0], yB[1], yB[2], (-M_PI+beta), TempVect2[0], TempVect2[1], TempVect2[2], nu, -PB[2]);
        }
        else         // Configuration II
        {   AngDisDispFSC_inDispHS(vA, yA[0], yA[1], yA[2], (      beta), TempVect2[0], TempVect2[1], TempVect2[2], nu, -PA[2]);
            AngDisDispFSC_inDispHS(vB, yB[0], yB[1], yB[2], (      beta), TempVect2[0], TempVect2[1], TempVect2[2], nu, -PB[2]);
        }
        // Calculate total Free Surface Correction to displacements in ADCS
        // Transform total Free Surface Correction to displacements from ADCS  to EFCS
        CoordTrans_inDispHS(d_vals, (vB[0]-vA[0]),(vB[1]-vA[1]),(vB[2]-vA[2]),Amat_t);

        DispVect[0] = d_vals[0];        DispVect[1] = d_vals[1];        DispVect[2] = d_vals[2];
    }
    
    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void AngDisDisp_inDispHS(double DispVect[3], double x, double y, double z, double alpha, double bx, double by, double bz, double nu)
{   // AngDisDisp calculates the "incomplete" displacements (without the  Burgers' function contribution) associated with an angular dislocation in an elastic full-space.
    double cosA,        sinA,       eta,        zeta,       r;
    double ux,          vx,         wx,         uy,         vy;
    double wy,          uz,         vz,         wz;
    //--------------------------------------
    cosA = cos(alpha);              sinA = sin(alpha);
    eta  = y*cosA-z*sinA;           zeta = y*sinA+z*cosA;
    r    = pow ((x*x +y*y +z*z),0.5);
    // Avoid complex results for the logarithmic terms
    if (zeta > r)           {   zeta = r;       }
    if (z    > r)           {   z    = r;       }
    //--------------------------------------
    ux = bx       /8.0/M_PI/(1.0-nu)*(x*y/r/(r-z)-x*eta/r/(r-zeta));
    vx = bx       /8.0/M_PI/(1.0-nu)*(eta*sinA/(r-zeta)-y*eta/r/(r-zeta)+y*y/r/(r-z)+(1.0-2.0*nu)*(cosA*log(r-zeta)-log(r-z)));
    wx = bx       /8.0/M_PI/(1.0-nu)*(eta*cosA/(r-zeta)-y/r-eta*z/r/(r-zeta)-(1.0-2.0*nu)*sinA*log(r-zeta));

    uy = by       /8.0/M_PI/(1.0-nu)*(x*x*cosA/r/(r-zeta)-                x*x/r/(r-z)-(1.0-2.0*nu)*(cosA*log(r-zeta)-log(r-z)));
    vy = by*x     /8.0/M_PI/(1.0-nu)*(  y*cosA/r/(r-zeta)-sinA*cosA/(r-zeta)-y/r/(r-z));
    wy = by*x     /8.0/M_PI/(1.0-nu)*(  z*cosA/r/(r-zeta)-cosA*cosA/(r-zeta)+1.0/r);

    uz = bz  *sinA/8.0/M_PI/(1.0-nu)*((1.0-2.0*nu)*log(r-zeta)-x*x/r/(r-zeta));
    vz = bz*x*sinA/8.0/M_PI/(1.0-nu)*(sinA/(r-zeta)-y/r/(r-zeta));
    wz = bz*x*sinA/8.0/M_PI/(1.0-nu)*(cosA/(r-zeta)-z/r/(r-zeta));

    DispVect[0] = ux+uy+uz;
    DispVect[1] = vx+vy+vz;
    DispVect[2] = wx+wy+wz;

    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void   AngDisDispFSC_inDispHS(double DispVect[3], double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a)
{   // AngDisDispFSC calculates the harmonic function contribution to the  displacements associated with an angular dislocation in an elastic  half-space.

    double sinB,        cosB,       cotB,       y3b,        z1b;
    double z3b,         r2b,        rb,         Fib,        v1cb1;
    double v2cb1,       v3cb1,      v1cb2,      v2cb2,      v3cb2;
    double v1cb3,       v2cb3,      v3cb3,      tempval1;


    sinB = sin(beta);               cosB = cos(beta);           cotB = cosB/sinB;
    y3b  = y3+2.0*a;                z1b  = y1*cosB+y3b*sinB;    z3b  = -y1*sinB+y3b*cosB;
    r2b  = y1*y1+y2*y2+y3b*y3b;     rb   = pow(r2b,0.5);

    tempval1 = (-(rb+y3b)* (cos(beta/2.0)/sin(beta/2.0)) +y1);                 
    Fib      = 2.0* atan(-y2/tempval1); //The Burgers' function

    v1cb1 = b1/4.0/M_PI/(1.0-nu)*(-2.0*(1.0-nu)*(1.0-2.0*nu)*Fib*pow(cotB,2.0)+(1.0-2.0*nu)*y2/(rb+y3b)*((1.0-2.0*nu-a/rb)*cotB-y1/(rb+y3b)*(nu+a/rb))+(1.0-2.0*nu)*y2*cosB*cotB/(rb+z3b)*(cosB+a/rb)+a*y2*(y3b-a)*cotB/pow(rb,3.0)+y2*(y3b-a)/(rb*(rb+y3b))*(-(1.0-2.0*nu)*cotB+y1/(rb+y3b)*(2.0*nu+a/rb)+a*y1/pow(rb,2.0))+y2*(y3b-a)/(rb*(rb+z3b))*(cosB/(rb+z3b)*((rb*cosB+y3b)*((1.0-2.0*nu)*cosB-a/rb)*cotB+2*(1.0-nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/pow(rb,2.0)));

    v2cb1 = b1/4.0/M_PI/(1.0-nu)*((1.0-2.0*nu)*((2.0*(1.0-nu)*pow(cotB,2.0)-nu)*log(rb+y3b)-(2.0*(1.0-nu)*pow(cotB,2.0)+1.0-2.0*nu)*cosB*log(rb+z3b))-(1.0-2.0*nu)/(rb+y3b)*(y1*cotB*(1.0-2.0*nu-a/rb)+nu*y3b-a+pow(y2,2.0)/(rb+y3b)*(nu+a/rb))-(1.0-2.0*nu)*z1b*cotB/(rb+z3b)*(cosB+a/rb)-a*y1*(y3b-a)*cotB/pow(rb,3.0)+(y3b-a)/(rb+y3b)*(-2.0*nu+1.0/rb*((1.0-2.0*nu)*y1*cotB-a)+pow(y2,2.0)/(rb*(rb+y3b))*(2.0*nu+a/rb)+a*pow(y2,2.0)/pow(rb,3.0))+(y3b-a)/(rb+z3b)*(pow(cosB,2.0)-1.0/rb*((1.0-2.0*nu)*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/pow(rb,3.0)-1.0/(rb*(rb+z3b))*(pow(y2,2.0)*pow(cosB,2.0)-a*z1b*cotB/rb*(rb*cosB+y3b))));

    v3cb1 = b1/4.0/M_PI/(1.0-nu)*(2.0*(1.0-nu)*(((1.0-2.0*nu)*Fib*cotB)+(y2/(rb+y3b)*(2.0*nu+a/rb))-(y2*cosB/(rb+z3b)*(cosB+a/rb)))+y2*(y3b-a)/rb*(2.0*nu/(rb+y3b)+a/pow(rb,2.0))+y2*(y3b-a)*cosB/(rb*(rb+z3b))*(1.0-2.0*nu-(rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)-a*y3b/pow(rb,2.0)));

    v1cb2 = b2/4.0/M_PI/(1.0-nu)*((1.0-2.0*nu)*((2.0*(1.0-nu)*pow(cotB,2.0)+nu)*log(rb+y3b)-(2.0*(1.0-nu)*pow(cotB,2.0)+1.0)*cosB*log(rb+z3b))+(1.0-2.0*nu)/(rb+y3b)*(-(1.0-2.0*nu)*y1*cotB+nu*y3b-a+a*y1*cotB/rb+pow(y1,2.0)/(rb+y3b)*(nu+a/rb))-(1.0-2.0*nu)*cotB/(rb+z3b)*(z1b*cosB-a*(rb*sinB-y1)/(rb*cosB))-a*y1*(y3b-a)*cotB/pow(rb,3.0)+(y3b-a)/(rb+y3b)*(2.0*nu+1.0/rb*((1.0-2.0*nu)*y1*cotB+a)-pow(y1,2.0)/(rb*(rb+y3b))*(2.0*nu+a/rb)-a*pow(y1,2.0)/pow(rb,3.0))+(y3b-a)*cotB/(rb+z3b)*(-cosB*sinB+a*y1*y3b/(pow(rb,3.0)*cosB)+(rb*sinB-y1)/rb*(2.0*(1.0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b)*(1.0+a/(rb*cosB)))));
                    
    v2cb2 = b2/4.0/M_PI/(1.0-nu)*(2.0*(1.0-nu)*(1.0-2.0*nu)*Fib*pow(cotB,2.0)+(1.0-2.0*nu)*y2/(rb+y3b)*(-(1.0-2.0*nu-a/rb)*cotB+y1/(rb+y3b)*(nu+a/rb))-(1.0-2.0*nu)*y2*cotB/(rb+z3b)*(1.0+a/(rb*cosB))-a*y2*(y3b-a)*cotB/pow(rb,3.0)+y2*(y3b-a)/(rb*(rb+y3b))*((1.0-2.0*nu)*cotB-2.0*nu*y1/(rb+y3b)-a*y1/rb*(1.0/rb+1.0/(rb+y3b)))+y2*(y3b-a)*cotB/(rb*(rb+z3b))*(-2.0*(1.0-nu)*cosB+(rb*cosB+y3b)/(rb+z3b)*(1.0+a/(rb*cosB))+a*y3b/(pow(rb,2.0)*cosB)));
                    
    v3cb2 = b2/4.0/M_PI/(1.0-nu)*(-2.0*(1.0-nu)*(1.0-2.0*nu)*cotB*(log(rb+y3b)-cosB*log(rb+z3b))-2.0*(1.0-nu)*y1/(rb+y3b)*(2.0*nu+a/rb)+2.0*(1.0-nu)*z1b/(rb+z3b)*(cosB+a/rb)+(y3b-a)/rb*((1.0-2.0*nu)*cotB-2.0*nu*y1/(rb+y3b)-a*y1/pow(rb,2.0))-(y3b-a)/(rb+z3b)*(cosB*sinB+(rb*cosB+y3b)*cotB/rb*(2.0*(1.0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b))+a/rb*(sinB-y3b*z1b/pow(rb,2.0)-z1b*(rb*cosB+y3b)/(rb*(rb+z3b)))));

    v1cb3 = b3/4.0/M_PI/(1.0-nu)*((1.0-2.0*nu)*(y2/(rb+y3b)*(1.0+a/rb)-y2*cosB/(rb+z3b)*(cosB+a/rb))-y2*(y3b-a)/rb*(a/pow(rb,2.0)+1.0/(rb+y3b))+y2*(y3b-a)*cosB/(rb*(rb+z3b))*((rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)+a*y3b/pow(rb,2.0)));
                    
    v2cb3 = b3/4.0/M_PI/(1.0-nu)*((1.0-2.0*nu)*(-sinB*log(rb+z3b)-y1/(rb+y3b)*(1.0+a/rb)+z1b/(rb+z3b)*(cosB+a/rb))+y1*(y3b-a)/rb*(a/pow(rb,2.0)+1.0/(rb+y3b))-(y3b-a)/(rb+z3b)*(sinB*(cosB-a/rb)+z1b/rb*(1.0+a*y3b/pow(rb,2.0))-1.0/(rb*(rb+z3b))*(pow(y2,2.0)*cosB*sinB-a*z1b/rb*(rb*cosB+y3b))));
                    
    v3cb3 = b3/4.0/M_PI/(1.0-nu)*(2.0*(1.0-nu)*Fib+2.0*(1.0-nu)*(y2*sinB/(rb+z3b)*(cosB+a/rb))+y2*(y3b-a)*sinB/(rb*(rb+z3b))*(1.0+(rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)+a*y3b/pow(rb,2.0)));

    DispVect[0] = v1cb1+v1cb2+v1cb3;
    DispVect[1] = v2cb1+v2cb2+v2cb3;
    DispVect[2] = v3cb1+v3cb2+v3cb3;

    return;
}

