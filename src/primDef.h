// code for primitive variable definition:
 
  double rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E;

  rho = Prim[0];
  v1 = Prim[1];
  v2 = Prim[2];
  theta  = Prim[3];
  T11 = Prim[4];
  T12 = Prim[5];
  T22 = Prim[6];
  q1 = Prim[7];
  q2 = Prim[8];
  vs = v1*v1+v2*v2;
  E = 0.5*rho*vs+1.5*rho*theta;
