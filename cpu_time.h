#ifndef CPU_TIME_H
#define CPU_TIME_H

double cpu_time ( ); 
 
#endif 

double cpu_time ( )

//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;
  
  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;
  
  return value;
}
