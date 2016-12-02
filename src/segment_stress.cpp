#include <apf.h>

namespace Dislocation {

static apf::Matrix3x3 outer(
    apf::Vector3 const& a,
    apf::Vector3 const& b) {
  apf::Matrix3x3 r;
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    r[i][j] = a[i]*b[j];
  return r;
}

static apf::Matrix3x3 triple_symmetric(
    apf::Vector3 const& a,
    apf::Vector3 const& b,
    apf::Vector3 const& c) {
  auto d = apf::cross(a,b);
  auto e = outer(d,c);
  auto f = outer(c,d);
  return (e+f)*0.5;
}

static apf::Matrix3x3 identity() {
  apf::Matrix3x3 I;
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    I[i][j] = 0.0;
  for (int i=0; i < 3; ++i)
    I[i][i] = 1.0;
  return I;
}

static apf::Matrix3x3 compute_point_stress(
    double mu,
    double nu,
    apf::Vector3 const& r,
    apf::Vector3 const& bp,
    apf::Vector3 const& tp,
    apf::Vector3 const& rp) {

  auto R = r-rp;
  auto Rm = R.getLength();
  auto Lp = R*tp;
  auto rho = R - tp*Lp;
  auto Y = tp*(Lp+Rm) + rho;
  auto Y2 = 2.0*Rm*(Lp+Rm);

  auto T1 = triple_symmetric(bp,Y,tp);
  auto T2 = triple_symmetric(bp,tp,Y);
  auto I = identity();

  auto s1 = mu/(apf::pi*Y2);
  auto s2 = 1.0/(1.0-nu);
  auto s3 = (bp*apf::cross(Y,tp))/(2.0*(1.0-nu));
  auto s4 = 2.0/Y2;
  auto s5 = Lp/Rm;

  apf::Matrix3x3 stress;
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    stress[i][j] =
      s1*(T1[i][j] - s2*T2[i][j] - s3*(I[i][j] + tp[i]*tp[j] +
            s4*(rho[i]*Y[j] + rho[j]*Y[i] + s5*Y[i]*Y[j])));
  return stress;
}

apf::Matrix3x3 compute_segment_stress(
    double mu,
    double nu,
    apf::Vector3 const& A,
    apf::Vector3 const& B,
    apf::Vector3 const& bp,
    apf::Vector3 const& r) {
  auto tp = B-A;
  tp = tp/tp.getLength();
  auto sB = compute_point_stress(mu, nu, r, bp, tp, B);
  auto sA = compute_point_stress(mu, nu, r, bp, tp, A);
  auto s = sB-sA;
  return s;
}

}
