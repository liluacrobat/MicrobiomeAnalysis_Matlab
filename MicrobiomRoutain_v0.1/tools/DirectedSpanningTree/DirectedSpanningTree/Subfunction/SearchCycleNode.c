# include "mex.h"
# include "math.h"
/*function TreeMatric = SearchCycleNode( ClusterNode,RootSharedNode,TreeMatric,OriginalCostMatric ) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int p,q,VarLen;
   double *Temp;
   VarLen = mxGetN( prhs[0] ) > mxGetM( prhs[0] ) ? mxGetN( prhs[0] ):mxGetM( prhs[0] );
   int ClusterNode[ VarLen ];
   
   /******************ClusterNode**********************/
   /*printf("\nVarLen is %d \nCluster nodes are:\n",VarLen);*/   
   Temp = mxGetPr( prhs[0] );
   for( p=0;p<VarLen;p++ ){
        ClusterNode[p] = (int)Temp[p]-1;   
        /*printf("%d ",ClusterNode[p]);*/
   }
   /*printf("\n");*/
   
   /********************RootSharedNode*********************/
   int RootLen = mxGetM( prhs[1] ) > mxGetN( prhs[1] ) ? mxGetM( prhs[1] ) : mxGetN(prhs[1]);
   /*printf("\nThe RootLen is %d\n",RootLen );*/
   int RootSharedNode[RootLen];
   Temp = mxGetPr(prhs[1]);
   /*printf("\nThe RootSharedNodes are\n",RootLen );*/
   for( p=0;p<RootLen;p++ ){
        RootSharedNode[p] = (int)Temp[p]-1;
        /*printf("%d ", RootSharedNode[p]+1 );*/
   }
  
   /****************TreeMatric************************/
   int Dim = mxGetM(prhs[2]);
   double **TreeMatric = mxCalloc( Dim, sizeof(double *) );
   for( p=0;p<Dim;p++ )
       TreeMatric[ p ] = mxCalloc( Dim, sizeof(double) );
   
   /* printf("\nTree matric is \n"); */
   Temp = mxGetPr( prhs[2] );
   for(p=0;p<Dim;p++){
       for(q=0;q<Dim;q++){
           TreeMatric[p][q] = Temp[ p+q*Dim ];
          /* printf("%.2f ",TreeMatric[p][q]);*/
       }
       /* printf("\n");*/
   }
 
   /*********************OriginalCostMatric***************/
   double **OriginalCostMatric = mxCalloc( Dim, sizeof(double *) );
   for( p=0;p<Dim;p++ )
       OriginalCostMatric[ p ] = mxCalloc( Dim, sizeof(double) );
   
   /* printf("\nOriginalcostMatric is\n");*/
   Temp = mxGetPr(prhs[3]);
   for(p=0;p<Dim;p++){
       for(q=0;q<Dim;q++){
           OriginalCostMatric[p][q] = Temp[ p+q*Dim ];
          /* printf("%.2f ",OriginalCostMatric[p][q]);*/
       }
       /* printf("\n");*/
   }
 
   
   /****DFS finds start point of the loop, whose nodes saved in ClusterNode****/
   int i=0,Run=1,Start,End,Len,Descendant;
   int Visited[Dim],List[Dim];
   while( i<VarLen && Run ==1 ){
          Start = ClusterNode[i];
          i++; p = 0;
          while( p<Dim ){ /*look for the parent of the node Start */
              if( TreeMatric[p][Start]!=0.0 ) break;
              p++;
          }          
          if( p<Dim ){
              End = p;
              for(p=0;p<Dim;p++){ Visited[p]=0; List[p]=0; }
              Visited[ Start ] = 1; Len=0;
              for( p=0;p<Dim;p++ )
                   if( TreeMatric[Start][p] != 0.0){
                       List[Len]=p; Len++;          
                   }
              while( Len>0 && Run==1 ){
                    Len--;
                    Descendant = List[ Len ];
                    if( Descendant==End ){  /*find the cycle node, if one of the descendant meets the node End */
                        Run = 0; 
                        break;
                    }
                    if( Visited[ Descendant ] == 0 ){
                        Visited[ Descendant ] = 1;
                        for( p=0;p<Dim;p++ )
                          if( TreeMatric[Descendant][p] != 0.0 ){
                              List[Len]=p; Len++;
                          }
                     }
              }
         }
   }
   /* printf("\nStart node is %d\n",Start+1);*/
   
   /******************nodes in cycle are saved in CycleNode**************/
   int CycleNode[ VarLen ];
   double SumWeight = 0.0;
   CycleNode[0] = Start; Len = 1;
   while( 1 ){
       for( p=0;p<Dim;p++ )
            if( TreeMatric[ p ][ CycleNode[Len-1] ] != 0.0 ) break;
       SumWeight += OriginalCostMatric[ p ][ CycleNode[Len-1] ];
       if( p==Start ) break;
       CycleNode[ Len ] = p; 
       Len++;
   }
   
   for( p=1;p<(int)ceil(Len/2)+1;p++ ){
        q = CycleNode[ p ];
        CycleNode[ p ] = CycleNode[ Len-p ];
        CycleNode[ Len-p ] = q;
   }   
   /*
   printf("\nLen is %d \nSumWeight is %f\n CycleNode is \n",Len,SumWeight);
   for( p=0;p<Len;p++ )
        printf( "%d ",CycleNode[ p ]+1);
   printf("\n");   
   */
   
  /******************put detected cycle weights in CycleWeight**************/
   double CycleWeight[Len];
   CycleWeight[0] = SumWeight - OriginalCostMatric[ CycleNode[Len-1] ][ CycleNode[0] ];
   for( p=1;p<Len;p++ )
        CycleWeight[p] = SumWeight - OriginalCostMatric[ CycleNode[p-1] ][ CycleNode[p] ];
/*
   printf("CycleWeight is\n");
   for( p=0;p<Len;p++ )
        printf("%f ",CycleWeight[p] );
*/
   
   /*******************look for the optimal node for cutting**************/
   int SearchRootCluster, CurNode = -1;
   double LocalMax = -mxGetInf();   
   /*printf("RootLen and VarLen are (%d,%d) \n",RootLen,VarLen );*/
   for( p=0;p<RootLen;p++ )
       for( q=0;q<VarLen;q++ ){
           if( OriginalCostMatric[ RootSharedNode[p] ][ ClusterNode[q] ] != 0.0 ){
               /*printf("Test arc (%d,%d) \n",RootSharedNode[ p ]+1,ClusterNode[ q ]+1 );*/
               *Temp = OriginalCostMatric[ RootSharedNode[p] ][ ClusterNode[q] ]; 
               CurNode = ClusterNode[q];
               while( 1 ){
                   for( i=0;i<Dim;i++ )
                        if( TreeMatric[ i ][ CurNode ] != 0.0 ) break;                
                   *Temp += OriginalCostMatric[ CurNode ][ i ] - TreeMatric[ i ][ CurNode ];
                   /*Move the the parent node*/
                   CurNode = i;
                   /*Test whether current node is inside the cycle or not*/
                   i = 0;
                   while( i < Len ){
                       if( CycleNode[i] == CurNode ) break;
                       i++;
                   }                   
                   /*If the current node is inside the cycle, sum up the weight */
                   if( i < Len ){
                       *Temp += CycleWeight[i];
                       break;
                   }
               }
               /* printf("By adding (%d,%d), the weight is %f \n",RootSharedNode[ p ]+1,ClusterNode[ q ]+1,*Temp );*/   
               /* find the maximal weight*/
               if( LocalMax < *Temp){
                   LocalMax = *Temp;
                   Start    = p;
                   End      = q;
               /* printf("Add the arc (%d,%d), with LocalMax: %f \n",RootSharedNode[ Start ]+1,ClusterNode[ End ]+1,LocalMax );*/
               }
           }
       }
   /* printf("\nDetermine to add the arc (%d,%d)\n",RootSharedNode[ Start ]+1,ClusterNode[ End ]+1 ); */
      
   /************** break the cycle and reconstruct the matrix *************/
   for(i=0;i<Dim;i++) Visited[ i ] = 0;
   CurNode = ClusterNode[End]; 
   while( 1 ){
         /* printf( "Current node is %d\n",CurNode+1 ); */
         i = 0;
         while( i<Dim ){
              if( TreeMatric[ i ][ CurNode ] != 0.0 && Visited[ i ] == 0 ) break;
              i++;
         }
        /* printf("break arc (%d,%d)\n",i+1,CurNode+1);  
         printf("value is %.1f\n",TreeMatric[ 8 ][ 7 ] );*/
         if(i==Dim)break;
         
         TreeMatric[ i ][ CurNode ] = 0.0;
         TreeMatric[ CurNode ][ i ] = OriginalCostMatric[ CurNode ][ i ];                 
         Visited[ CurNode ] = 1;
         CurNode = i;
   }
   
   TreeMatric[ RootSharedNode[ Start ]][ClusterNode[ End ] ] = OriginalCostMatric[ RootSharedNode[ Start ] ][ ClusterNode[ End ] ];   
   
   /**************** return the recontructed matrix back******************/
   plhs[0] = mxCreateDoubleMatrix( Dim,Dim,mxREAL );
   Temp = mxGetPr( plhs[0] );
   for( p=0;p<Dim;p++ )
        for( q=0;q<Dim;q++ ){
             Temp[ p+q*Dim ] = TreeMatric[ p ][ q ];
        }
   
}
   