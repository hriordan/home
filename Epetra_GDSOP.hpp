//@HEADER
//Epetra GDS Implementation, Basic Operations - Vector, CrsMatrix, Map, and 
// SerialDenseSolver Classes 
//Henry Riordan - LSSG Research and development intern, 2013
//@HEADER
#ifndef EPETRA_GDSOP_HPP
#define EPETRA_GDSOP_HPP

#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>

#ifdef HAVE_gds


#include <gds.h>
#include <mpi.h>

#endif



namespace Epetra{


#ifdef HAVE_gds



  GDS_status_t Epetra_GDS_init(){

    GDS_thread_support_t provd_support;
    GDS_status_t ret = GDS_STATUS_OK;

    ret = GDS_init(NULL, NULL, GDS_THREAD_SINGLE, &provd_support);
    return ret;

  }


  GDS_status_t Epetra_GDS_finalize(){

    GDS_status_t ret;
    ret = GDS_finalize();
    return ret; 

  }

  //============================================================

  /*Epetra_Vector GDS functions */

  GDS_status_t Epetra_Vector_gds_alloc(Epetra_Vector &vec, GDS_gds_t *gds_v){
  
    GDS_size_t cts[1];
    int ndim = 1;
    GDS_size_t min_chunk[]={0};

    cts[0] = vec.GlobalLength();

    GDS_status_t ret = GDS_STATUS_OK;
    
      
    ret = GDS_alloc(ndim, cts, min_chunk, GDS_DATA_DBL, GDS_PRIORITY_HIGH, GDS_COMM_WORLD,MPI_INFO_NULL, gds_v);
  
      
    return ret;

  }


  GDS_status_t Epetra_Vector_gds_put(Epetra_Vector &vec, GDS_gds_t gds_v){

    GDS_size_t lo[1];
    GDS_size_t hi[1];
    GDS_size_t ld[] = {0};

    GDS_status_t ret = GDS_STATUS_OK;

    int ProcHasData = vec.MyLength();

    if ( ProcHasData > 0){

      lo[0] = vec.Map().MinMyGID();
      hi[0] = vec.Map().MaxMyGID();

      double *view;
      vec.ExtractView(&view); 

   
      ret = GDS_put(view, ld, lo, hi, gds_v);
 
 
    }

      GDS_wait(gds_v); //Ensure synchronization 
      
    return ret;    
  }


  GDS_status_t Epetra_Vector_gds_get(Epetra_Vector &vec, GDS_gds_t gds_v){

    GDS_size_t lo[1];
    GDS_size_t hi[1];


    GDS_size_t ld[] = {0};
    GDS_status_t ret= GDS_STATUS_OK;

    
    int ProcHasData = vec.MyLength(); 

    if (ProcHasData > 0){

      lo[0] = vec.Map().MinMyGID();
      hi[0] = vec.Map().MaxMyGID();

      double * view;
      vec.ExtractView(&view);

      ret = GDS_get(view, ld, lo, hi, gds_v);
      
    }
    GDS_wait_local(gds_v); 

    return ret;

  }


  GDS_status_t Epetra_Vector_gds_version_inc(GDS_gds_t gds_v){

    GDS_status_t ret;

    ret = GDS_version_inc(gds_v, 1, NULL, 0);
 
    return ret;
 
  }



  GDS_status_t Epetra_Vector_gds_version_dec(GDS_gds_t gds_v){

    GDS_status_t ret;

    ret = GDS_version_dec(gds_v, 1);

    return ret;
    
  }



  //========================================================

  /* Epetra_CrsMatrix GDS functions */
    
  GDS_status_t Epetra_CM_gds_alloc (Epetra_CrsMatrix &M, GDS_gds_t *gds_data,  GDS_gds_t * gds_col_ind, GDS_gds_t * gds_row_ptr){
    
    
    if (!(M.Filled()))
      {
	std::cout<<"The matrix is not fill complete yet!"<<std::endl;
	return GDS_STATUS_INVALID;
      }


  
    GDS_size_t cts_data[1], cts_col_ind[1], cts_row_ptr[1];
    GDS_size_t min_chunk[] = {0};
    int ndim = 1;

    cts_data[0] = M.NumGlobalNonzeros();

    GDS_status_t ret = GDS_STATUS_OK;

    ret = GDS_alloc(ndim, cts_data, min_chunk, GDS_DATA_DBL, GDS_PRIORITY_HIGH, GDS_COMM_WORLD, MPI_INFO_NULL, gds_data);                                                           
                                                                
    cts_col_ind[0] = cts_data[0];                                                          
    cts_row_ptr[0] = M.NumGlobalRows() + 1;                                                  
                                                                                           
    ret = GDS_alloc(ndim, cts_col_ind, min_chunk, GDS_DATA_INT, GDS_PRIORITY_HIGH, GDS_COMM_WORLD, MPI_INFO_NULL, gds_col_ind);                                                     
    ret = GDS_alloc(ndim, cts_row_ptr, min_chunk, GDS_DATA_INT, GDS_PRIORITY_HIGH, GDS_COMM_WORLD, MPI_INFO_NULL, gds_row_ptr);                                                     
                                                                               
    return ret; 

  }


  GDS_status_t Epetra_CM_gds_put(Epetra_CrsMatrix &M,  GDS_gds_t gds_data,  GDS_gds_t gds_col_ind, GDS_gds_t gds_row_ptr){

    if(!(M.Filled())){
      printf("Matrix not fill complete.\n");
      return GDS_STATUS_INVALID;
    }

    GDS_status_t ret = GDS_STATUS_OK;
                                                         
    int ProcHasData = M.Map().NumMyElements();

    GDS_size_t buff[1], offBuff[1];
    GDS_size_t lo_d[1], hi_d[1], lo_i[1], hi_i[1], lo_r[1], hi_r[1];
    GDS_size_t ld[] = {0};

    unsigned int size = M.NumMyNonzeros();
    unsigned int offsetSize = M.NumMyRows();


   
    /*set pointers*/
    double *vals; 
    int *indices; 
    int *indexOffset;  
    
    M.ExtractCrsDataPointers(indexOffset, indices, vals);

    buff[0] = size;
    
    MPI_Scan(buff, hi_d, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    lo_d[0] = hi_d[0] - size; 
    hi_d[0] -= 1; 
  
    if (ProcHasData > 0){
      ret = GDS_put(vals, ld, lo_d, hi_d, gds_data);
      GDS_wait(gds_data);
    }

    if (ret == GDS_STATUS_OK){    
    
      buff[0] = size;

      MPI_Scan(buff, hi_i, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      lo_i[0] = hi_i[0] - size;
      hi_i[0] -= 1;
 
      if(ProcHasData > 0)
	ret = GDS_put(indices, ld, lo_i, hi_i, gds_col_ind);
      GDS_wait(gds_col_ind);

    }

 
    if (ret == GDS_STATUS_OK){
   
      offBuff[0] = offsetSize; 
    
      MPI_Scan(offBuff, hi_r, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
   
      lo_r[0] = hi_r[0] - offsetSize;
      hi_r[0] -= 1; 
    
      if (ProcHasData > 0)
	ret = GDS_put(indexOffset, ld, lo_r, hi_r, gds_row_ptr);  
      GDS_wait(gds_row_ptr);
    
    }

    return ret;

  } 

  GDS_status_t Epetra_CM_gds_get(Epetra_CrsMatrix &M, GDS_gds_t gds_data, GDS_gds_t gds_col_ind, GDS_gds_t gds_row_ptr){

  

    if(!(M.Filled())){
      printf("Matrix not fill complete.\n");
      return GDS_STATUS_INVALID;
    }
    
    GDS_status_t ret = GDS_STATUS_OK;
    
    int ProcHasData = M.Map().NumMyElements();
    GDS_size_t buff[1];
    GDS_size_t lo_d[1], hi_d[1], lo_i[1], hi_i[1], lo_r[1], hi_r[1];
    GDS_size_t ld[] = {0};
 
    int size = M.NumMyNonzeros();
    int offsetSize = M.NumMyRows();
  
    double *vals;
    int *indices; 
    int *indexOffset; 
 
    M.ExtractCrsDataPointers(indexOffset, indices, vals); /*Equiv to a view of the data */
   
    buff[0] = size;

    MPI_Scan(buff, hi_d, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
 
    lo_d[0] = hi_d[0] - size;
    hi_d[0] -= 1; 

    if (ProcHasData > 0)
      ret = GDS_get(vals, ld, lo_d, hi_d, gds_data); 
    GDS_wait_local(gds_data);

    if (ret == GDS_STATUS_OK){
   
      buff[0] = size;

      MPI_Scan(buff, hi_i, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      lo_i[0] = hi_i[0] - size;
      hi_i[0] -= 1; 

      if (ProcHasData > 0)
	ret = GDS_get(indices, ld, lo_i, hi_i, gds_col_ind);
      GDS_wait_local(gds_col_ind);

    }


    if (ret == GDS_STATUS_OK){

      buff[0] = offsetSize;
 
      MPI_Scan(buff, hi_r,1,MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      lo_r[0] = hi_r[0] - offsetSize;
      hi_r[0] -= 1;
      
      if (ProcHasData > 0)
	ret = GDS_get(indexOffset, ld, lo_r, hi_r, gds_row_ptr);
      GDS_wait_local(gds_row_ptr);
    } 
  
  

    return ret; 
 
  } 

  GDS_status_t Epetra_CM_gds_version_inc(GDS_gds_t  gds_data, GDS_gds_t gds_col_ind, GDS_gds_t  gds_row_ptr ){

    GDS_status_t ret; 

    ret = GDS_version_inc(gds_data, 1, NULL, 0);
  
    if(ret != GDS_STATUS_OK)
      return ret; 

    ret = GDS_version_inc(gds_row_ptr, 1, NULL, 0);

    if(ret != GDS_STATUS_OK)
      return ret; 

    ret = GDS_version_inc(gds_col_ind, 1, NULL, 0);

    if(ret != GDS_STATUS_OK)
      return ret; 



    return GDS_STATUS_OK;

  }

  GDS_status_t Epetra_CM_gds_version_dec(GDS_gds_t  gds_data, GDS_gds_t gds_col_ind, GDS_gds_t  gds_row_ptr ){

    GDS_status_t ret; 

    ret = GDS_version_dec(gds_data, 1);
  
    if(ret != GDS_STATUS_OK)
      return ret; 

    ret = GDS_version_dec(gds_row_ptr, 1);

    if(ret != GDS_STATUS_OK)
      return ret; 

    ret = GDS_version_dec(gds_col_ind, 1);

    if(ret != GDS_STATUS_OK)
      return ret; 

    return GDS_STATUS_OK;

  }


  //===========================================================================

 //* helper gds functions for SerialDenseMatrices within SerialDenseSolver *//

  GDS_status_t Epetra_SDM_gds_alloc(Epetra_SerialDenseMatrix &M, GDS_gds_t *gds_v){
    
    int ndim = 1;
    GDS_size_t cts[1];
    GDS_size_t min_chunk[] = {0};
    GDS_status_t ret; 

    cts[0] = M.M() * M.N();

    ret = GDS_alloc(ndim, cts, min_chunk, GDS_DATA_DBL, GDS_PRIORITY_HIGH, GDS_COMM_WORLD, MPI_INFO_NULL, gds_v);

    return ret; 
  }

    
  GDS_status_t Epetra_SDM_gds_put(Epetra_SerialDenseMatrix &M, GDS_gds_t gds_v){
    
    GDS_size_t lo[1];
    GDS_size_t hi[1];
    GDS_size_t ld[] = {0};

   
   
    GDS_status_t ret = GDS_STATUS_OK;

    hi[0] = (M.M() * M.N()) - 1;
    lo[0] = 0; 
    
    double *x = M.A();

    ret = GDS_put(x, ld, lo, hi, gds_v);
    GDS_wait(gds_v);

    return ret; 
    
  }


  GDS_status_t Epetra_SDM_gds_get(Epetra_SerialDenseMatrix &M, GDS_gds_t gds_v){
  

    GDS_size_t lo[1];
    GDS_size_t hi[1];
    GDS_size_t ld[] = {0};

    GDS_status_t ret = GDS_STATUS_OK;

    hi[0] = (M.M() * M.N()) - 1;
    lo[0] = 0;

    double *x = M.A();

    ret = GDS_get(x, ld, lo, hi, gds_v);
    GDS_wait_local(gds_v);

    return ret;


  }

  GDS_status_t Epetra_SDM_gds_version_inc(GDS_gds_t gds_v){

    GDS_status_t ret;

    ret = GDS_version_inc(gds_v, 1, NULL, 0);

    return ret;
 
  }



  GDS_status_t Epetra_SDM_gds_version_dec(GDS_gds_t gds_v){

    GDS_status_t ret;

    ret = GDS_version_dec(gds_v, 1);

    return ret;
    
  }



 //============================================================
  /* SerialDenseSolver GDS wrap functions */

  /*require allocation for matrix, factored matrix, LHS and RHS; just to be safe for the moment */

  //The order of gds arguments is important, otherwise data will be mixed up.
  GDS_status_t Epetra_SDSolver_gds_alloc(Epetra_SerialDenseSolver &solver, GDS_gds_t *gds_matrix, GDS_gds_t *gds_fmatrix, GDS_gds_t *gds_lhs, GDS_gds_t *gds_rhs){

    GDS_status_t ret;
    
    ret = Epetra_SDM_gds_alloc(*solver.Matrix(), gds_matrix);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_alloc(*solver.FactoredMatrix(), gds_fmatrix);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_alloc(*solver.LHS(), gds_lhs);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_alloc(*solver.RHS(), gds_rhs);
 
    return ret;

  }

  
  GDS_status_t Epetra_SDSolver_gds_put(Epetra_SerialDenseSolver &solver, GDS_gds_t gds_matrix, GDS_gds_t gds_fmatrix, GDS_gds_t gds_lhs, GDS_gds_t gds_rhs){

    GDS_status_t ret;

    ret = Epetra_SDM_gds_put(*solver.Matrix(),gds_matrix);
    
    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_put(*solver.FactoredMatrix(),gds_fmatrix);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_put(*solver.LHS(),gds_lhs);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_put(*solver.RHS(),gds_rhs);

    return ret; 

  }


  GDS_status_t Epetra_SDSolver_gds_get(Epetra_SerialDenseSolver &solver, GDS_gds_t gds_matrix, GDS_gds_t gds_fmatrix, GDS_gds_t gds_lhs, GDS_gds_t gds_rhs){

    GDS_status_t ret; 

    ret =  Epetra_SDM_gds_get(*solver.Matrix(), gds_matrix);

    if (ret == GDS_STATUS_OK)
      ret =  Epetra_SDM_gds_get(*solver.FactoredMatrix(), gds_fmatrix);
  
    if (ret == GDS_STATUS_OK)    
      ret =  Epetra_SDM_gds_get(*solver.RHS(), gds_rhs);
      
    if (ret == GDS_STATUS_OK)
      ret =  Epetra_SDM_gds_get(*solver.LHS(), gds_lhs);

    return ret; 

  }

  GDS_status_t Epetra_SDSolver_gds_v_inc(GDS_gds_t gds_matrix, GDS_gds_t gds_fmatrix, GDS_gds_t gds_lhs, GDS_gds_t gds_rhs){

    GDS_status_t ret; 

    ret = Epetra_SDM_gds_version_inc(gds_matrix);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_inc(gds_fmatrix);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_inc(gds_lhs);
    
    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_inc(gds_rhs);

    return ret; 

  }

  GDS_status_t Epetra_SDSolver_gds_v_dec(GDS_gds_t gds_matrix, GDS_gds_t gds_fmatrix, GDS_gds_t gds_lhs, GDS_gds_t gds_rhs){

    GDS_status_t ret; 

    ret = Epetra_SDM_gds_version_dec(gds_matrix);
   
    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_dec(gds_fmatrix);
    
    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_dec(gds_lhs);

    if (ret == GDS_STATUS_OK)
      ret = Epetra_SDM_gds_version_dec(gds_rhs);

    return ret; 

  }
 
  //==============================================

  /* Epetra_Map GDS functionality */

  /*Core concept of this implementation, unique to map due to the class's data access restrictions, is the reconstruction of the map entirely in GDS getting */
  
  GDS_status_t Epetra_Map_gds_alloc(Epetra_Map &Map, GDS_gds_t *gds_elem/* ,  GDS_gds_t *gds_point, GDS_gds_t *gds_size, GDS_gds_t *gds_pointTo */){
    
    int ndim = 1;
    GDS_size_t min_chunk[] = {0};
    GDS_size_t cts[1];


    GDS_status_t ret = GDS_STATUS_OK;

    cts[0] = Map.NumGlobalElements();

    
    ret = GDS_alloc(ndim, cts, min_chunk, GDS_DATA_INT, GDS_PRIORITY_HIGH, GDS_COMM_WORLD, MPI_INFO_NULL, gds_elem);
   
    return ret; 
 
  }





  GDS_status_t Epetra_Map_gds_put(Epetra_Map &Map, GDS_gds_t gds_elem){


    GDS_status_t ret = GDS_STATUS_OK;
    GDS_size_t hi[1], lo[1];
    GDS_size_t ld[] = {0};

    int *elems = Map.MyGlobalElements();
    GDS_size_t buff[1];
    int size = Map.NumMyElements();
    
    buff[0] = size; 

    MPI_Scan(buff, hi, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    lo[0] = hi[0] - size;
    hi[0] -= 1;
    
    ret = GDS_put(elems, ld, lo, hi, gds_elem);
    GDS_wait(gds_elem);

    return ret;

  }

  /*
  IN: NumGlobalElements. Needed in map construction
  IN: NumMYElements. Ditto.
  IN: gds object
  IN: Comm object.
  OUT: Map! 
  */
  //TODO: use map object to extract eleme num details? But is the Map to be trusted in this context?
 
  GDS_status_t Epetra_Map_gds_get( Epetra_Map &Map, int NumGlobalElements, int NumMyElements ,GDS_gds_t gds_v, Epetra_MpiComm Comm){

    GDS_status_t ret = GDS_STATUS_OK;
    GDS_size_t hi[1], lo[1];
    GDS_size_t ld[] = {0};

    int elBuff[NumMyElements];
    int size = NumMyElements; 
    GDS_size_t buff[1]; 

    buff[0] = size; 
    
    MPI_Scan(buff, hi, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    
    lo[0] = hi[0] - size; 
    hi[0] -= 1; 
    
    ret = GDS_get(elBuff, ld, lo, hi, gds_v);   //elBuff now has globalElements needed for reconstruction
    GDS_wait_local(gds_v);

    
    Epetra_Map New(NumGlobalElements, NumMyElements, elBuff, 0, Comm);

    Map = New;

    return ret; 

  }

  /*IDEA: have 1-length gds's for each of the constructor params? So you don't need to preserve unique map param variables in main code?
    i.e. gds_NGE, gds_NME, etc. 
    * might be a waste of computation, though.
  */


  /* Class-particular versioning functions unnecessary in this implementation */ 


  //===============================================


#endif

  }

#endif
