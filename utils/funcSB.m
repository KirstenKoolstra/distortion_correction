% This function calculates b = Ax
function vec_out=funcSB(vec_in)

global RcEc ER mu lambda

temp = ER*vec_in;
temp2 = RcEc*temp;
vec_out = mu*temp2 + lambda*( Dxt(Dx(vec_in)) + Dyt(Dy(vec_in))) ;

end