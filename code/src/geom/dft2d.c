// This is the DFT2D function, it's really just the same as the 1D versions of DFT which was explained in an earlier section, only with an extra loop to do the calculation for every row or column, and repeated twice (for the columns, and for the rows), and this separately for every color component. The function take now two parameters for the dimensions: m and n, because an image is 2D. The parameter inverse can be set to true to do the inverse calculations instead. Again *gRe and *gIm are the input arrays, and *GRe and *GIm the output arrays. In this 2D version, the values are divided through m for the normal DFT and through n for the inverse DFT. There are also definitions where you should divide through m*n in one direction and through 1 in the other, or though sqrt(m*n) in both, but it doesn't really matter which one you pick, as long as the forward FT and the inverse together give the original image back. 
// source: http://student.kuleuven.be/~m0216922/CG/fourier.html#dft2d

#include "dft2d.h"

void dft2d(int n, int m, int inverse, double *gRe, double *gIm, double *GRe, double *GIm)
{
  double pi = 3.1415926535897932384626433832795;
  int nns, w, x, y;
  double *Gr2, *Gi2;
  double a, ca, sa; 
  double fm,fn;
  
  nns=n*m;
  //printf("n=%d, m=%d, nns=%d\n", n, m, nns);
  Gr2 = dvector(1,nns);   //temporary buffers
  Gi2 = dvector(1,nns);

  //calculate the fourier transform of the columns
  for(x = 0; x < n; x++)
  {
//    print(" % done",8, 0, RGB_White, 1);
//    print(50 * x / n, 0, 0, RGB_White, 1);
//    if(done()) end();
//    redraw();
//    //This is the 1D DFT:
    for(w = 0; w < m; w++)
    {
      Gr2[m * x + w+1 ] = Gi2[m * x + w+1 ] = 0;
      for(y = 0; y < m; y++)
      { fm=m;
        a= 2 * pi * w * y /fm;
        if(inverse)a = -a;
        ca = cos(a);
        sa = sin(a);
        Gr2[m * x + w+1] += gRe[m * x + y+1] * ca - gIm[m * x + y+1 ] * sa;
        Gi2[m * x + w+1] += gRe[m * x + y+1] * sa + gIm[m * x + y+1 ] * ca;    
      }
   
    }
  }
  //calculate the fourier transform of the rows
  for(y = 0; y <m; y++)
  {
//    print(" % done",8, 0, RGB_White, 1);
//    print(50 + 50 * y / m, 0, 0, RGB_White, 1);
//    if(done()) end();
//    redraw();
//    //This is the 1D DFT:
    for(w = 0; w < n; w++)
    {
      GRe[m * w + y+1 ]= GIm[m * w + y+1 ] = 0;
      for(x = 0; x < n; x++)
      { fn=n;
        a = 2 * pi * w * x / fn;
        if(inverse)a = -a;
        ca = cos(a);
        sa = sin(a);
        GRe[m * w + y+1] += Gr2[m * x + y+1] * ca - Gi2[m * x + y+1] * sa;
        GIm[m * w + y+1] += Gr2[m * x + y+1] * sa + Gi2[m * x + y+1] * ca;
      }
      if(inverse)
      {
        GRe[m * w + y+1] /= n;
        GIm[m * w + y+1] /= n;
      }
      else
      {
        GRe[m * w + y+1] /= m;
        GIm[m * w + y+1] /= m;
      }
    }
  }

  free_dvector(Gr2,1,nns);   //temporary buffers
  free_dvector(Gi2,1,nns);
}
