//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/examples/swe/sweFort.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: F77 external declarations for SAMRAI shallow water eqns ex.
//

extern "C" {

  void initflow_(
  const int& , 
  const double*, const double*, const double*,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  double*, double*, 
  double*, double*);

  void setbathy_(
  const int& , 
  const double*, const double*, const double*,
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  double*);

  void tagwd_(
  const double*, 
  const int& , const int& ,
  const int& , const int& ,
  const int& , const int& ,
  double*, double*, double*,
  double*, double*);
  

  void   calcdt_(const double*,
  const int& , const int& , 
  const int& , const int& ,
  const int& , const int& ,
  const double*, 
  const double*,
  const double*, 
  double&);

  void george_flux_(
  const int&, const double&, const double*,
  const int&, const int&, const int&, const int&,
  double*, double*, double*, 
  double*, double*, double*, double*);

  void consdiff_(
  const double*, 
  const int&, const int&, const int&, const int&,
  double*, double*, double*, double*,
  double*, double*, double* );

  void flux_sed_(
  const int&, const double&, const double*,
  const double*, const double*,
  const int&, const int&, const int&, const int&,
  double*, double*, double*, 
  double*, double* );

  void consdiff_sed_(
  const double*, 
  const int&, const int&, const int&, const int&,
  double*, double*, 
  double*, double* );

  void friction_(
  const double*, const double&,
  const int&, const int&, const int&, const int&,
  double*, double*, double*);

  void linfriction_(
  const double*, const double&,
  const int&, const int&, const int&, const int&,
  double*, double*, double*);

  void drycheck_(
  const double*, const double&,
  const int&, const int&, const int&, const int&,
  double*, double*);

  void junkprobe_(
  const int&, const double*, const double*, const double* ,const double&,
  const int&, const int&, const int&, const int&, 
  double*, double*, double*);

  void detectgrad_(
     const double&,
     const int& , const int& , 
     const int& , const int& , 
     const int& , const int& , 
     const int& , const int& ,
     const int& , const int& ,
     const double* , const double*, const double*,
     const double& , 
     const int&, const int&,
     const double*,
     int* , int* );

   void c2f_(const int& , const int& , const int&, const int& , const int&,
       const int&, const double&, const double&, const int&, const int&,
       const int&, const double&, const double&, const double& );

}
