inline double clamp_log(double x) {
  const double eps = 1e-10;
  return std::log(std::max(x, eps));
}