#include <vw/InterestPoint/InterestData.h>
#include <algorithm>

namespace vw {
namespace ip {

int cull_interest_points(std::vector<InterestPoint> &points, int num_points) {
	std::sort(points.begin(), points.end());
	if ((num_points >= 0) && (num_points < points.size()))
		points.erase(points.begin() + num_points, points.end());
	return points.size();
}

}} //vw::ip
