
const double eps = 1e-8;
const double PI = acos(-1);

bool zero(double d = 0){
	return abs(d) < eps;
}
struct point{
	double x, y;
	// Distance Squared
	double dist2(point p = {0, 0}){		
		return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
	}
	// Distance, by default from origin
	double dist(point p = {0, 0}) {
		return sqrt(dist2(p));
	}

	point operator + (const point & p) const {
		return {x + p.x, y + p.y};
	}

	point operator - (const point & p) const {
		return {x - p.x, y - p.y};
	}
	point operator * (const double mul) const {
		return {x * mul, y * mul};
	}

	bool operator < (const point& p) const{
		if(p.x == x) return y < p.y;
		return x < p.x;
	}
	bool operator == (const point & p) const{
		return x == p.x && y == p.y;
	}

	// return a point after rotating b by angle t wrt this
	point rotate(point b, double t){
		double xx = b.x - x;
		double yy = b.y - y;
		return {xx * cos(t) - yy * sin(t) + x,
				yy * cos(t) + xx * sin(t) + y};
	}
};

double cross(point a, point b){
	return a.x * b.y - a.y * b.x;
}

// f = 1 means 'signed'. Twice the area of triangle 
double area(point a, point b, point c, bool f = 0){
	double ar = a.x * b.y - b.x * a.y;
	ar += b.x * c.y - c.x * b.y;
	ar += c.x * a.y - a.x * c.y;
	if(!f) ar = abs(ar);
	return ar;
}

struct line{
	double a, b, c;
	// returns points of intersection of the lines. 
	// If they overlap then returns two origins
	vector<point> intersect(line l){
		double p = l.a, q = l.b, r = l.c;
		if(zero(b * p - a * q)){
			if(zero(b * r - c * q)) return {(point){0, 0}, (point){0, 0}};
			return {};
		}
		return {(point){(b*r-c*q)/(a*q-b*p),(c*p-a*r)/(a*q-b*p)}};
	}
	// fixed this: added abs for distance
	double dist(point p){
		double num = a * p.x + b * p.y + c;
		double den = sqrt(a * a + b * b);
		num = abs(num);
		return num / den;
	}

	// line perpendicular to this line passing through point p
	line perp(point p){
		return (line){b, -a, a * p.y - b * p.x};
	}
	// get two points at distance = d from a given point on the same line
	vector<point> points_at_d(point p, double d){
		point dir = {-b, a};
		return {p + dir * (d / dir.dist()) , p - dir * (d / dir.dist())};
	}

	// return 1 on same side otherwise 0
	// same side means non of the points should lie on the line
	bool same_side(point p1, point p2){
		double p = p1.x * a + b * p1.y + c;
		double q = p2.x * a + b * p2.y + c;
		if(p < 0) p *= -1, q *= -1;
		return q > 0;
	}
};

struct line_segment{
	point a, b;
	// returns perpendicular bisector line of the line segment
	line perp_bisec(){
		point mid = (a + b) * 0.5;
		double l = a.x - b.x;
		double m = a.y - b.y;
		double n = -(mid.y * m + l * mid.x);
		return (line){l, m, n};
	}
	// get standard line equation from two point form
	line get_line(){
		double l = b.y - a.y;
		double m = a.x - b.x;
		double n = -(a.y * m + a.x * l);
		return (line){l, m, n};
	}
	// Check if a point p lies inside, including end points, a line segment
	bool inside(point c){
		if(!zero((b.y - a.y) * (c.x - a.x) - (b.x - a.x) * (c.y - a.y))) return 0;
		if(c.x < min(a.x, b.x) || c.x > max(a.x, b.x)) return 0;
		if(c.y < min(a.y, b.y) || c.y > max(a.y, b.y)) return 0;
		return 1;
	}

	// returns true if line segments intersects or touches or overlaps.
	bool intersect(line_segment s){
		line l1 = get_line();
		line l2 = s.get_line();
		vector<point> v = l1.intersect(l2);
		if(v.empty()) return 0;
		if(v.size() > 1) {
			return (s.inside(a) || s.inside(b) || inside(s.a) || inside(s.b));
		}
		return s.inside(v[0]) && inside(v[0]);
	}
};

struct circle
{
	double xc, yc, r;
	// returns theta in [0, 2PI)
	double get_theta(point p){ //Note : p need not lie in the circle
		double t = atan2(p.y - yc, p.x - xc);
		if(t < 0) t += 2.0 * PI;
		return t;
	}
	// return a vector of points which intersect with the line
	vector<point> line_intersect(line l){
		double a = l.a, b = l.b, c = l.c;
		vector<point> ans;
		double d = l.dist({xc, yc});
		if(d > r + eps) return ans;
		// tangent
		line bis = l.perp((point){xc, yc});
		if(zero(d - r)) {
			return bis.intersect(l);
		}
		// secant
		point mid = bis.intersect(l)[0];
		double dd = r * r - d * d;
		dd = sqrt(dd);
		return l.points_at_d(mid, dd);
	}

	// return vector of points of intesection of the circles
	vector<point> intersect(circle c){
		double x1 = xc, y1 = yc, r1 = r;
		double x2 = c.xc, y2 = c.yc, r2 = c.r;

		if(r1 < r2) return c.intersect((circle){x1, y1, r1});

		point c1 = {x1, y1};
		point c2 = {x2, y2};

		double d = c1.dist(c2);

		// one intesecting point
		if(zero(c1.dist2(c2) - (r1+r2)*(r1+r2)) || zero(c1.dist2(c2) - (r1-r2)*(r1-r2))){
			return {{x1 + (x2-x1) * r1 / d,
					y1 + (y2-y1) * r1 / d}};
		}

		// no intesection point 
		if(c1.dist2(c2)+eps > (r1+r2)*(r1+r2) || c1.dist2(c2) < (r1-r2)*(r1-r2) + eps){
			return {};
		}

		double Cos_t = (r1*r1 + d*d - r2*r2)*0.5/d/r1;
		double t = acos(Cos_t);

		point p = {x1 + (x2-x1)*r1/d,
					y1 + (y2-y1)*r1/d};

		point p1 = c1.rotate(p, t);
		point p2 = c1.rotate(p, -t);
		return {p1, p2};
	}

	// return true if line segment l intesects / touches with the circle
	bool intersect(line_segment l){
		if(l.a.dist2({xc, yc}) < r*r - eps || l.b.dist2({xc, yc}) < r*r - eps) 
			return 1;
		return line_intersect(l.get_line()).size() > 1;
	}

	// return points where tangent from point p touch the circles
	// p may be anywhere in the plane. All cases handled
	vector <point> tangent_points(point p){
		if(zero(p.dist2({xc, yc}) - r*r)) return {p};
		if(p.dist2({xc, yc}) < r*r + eps) return {};
		double d = p.dist({xc, yc});
		point cen = {xc, yc};
		double t = acos(min((double)1, r / d));
		double x1 = xc + (p.x - xc) * r / d;
		double y1 = yc + (p.y - yc) * r / d;
		point pp = {x1, y1};
		vector<point> ans;
		ans.push_back(cen.rotate(pp, t));
		ans.push_back(cen.rotate(pp, -t));
		
		return ans;
	}

	// returns cross tangents between the circles as vector of line segments
	// All cases handled. Cross tangent is not possible if circles overlap.
	// In case of external-touch, line segment of zero length passing through point of contact is returned.
	vector<line_segment> cross_tangents(circle c){
		double x1 = xc, y1 = yc, r1 = r;
		double x2 = c.xc, y2 = c.yc, r2 = c.r;

		point c1 = {x1, y1};
		point c2 = {x2, y2};

		// one tangent : returns a zero length line segment through the point of contact 
		if(zero(c1.dist2(c2) - (r1+r2)*(r1+r2))){
			c1 = intersect(c)[0];
			return {{c1, c1}};
		}

		// circles intesecting / touching
		if(c1.dist2(c2) <= (r1+r2)*(r1+r2) + eps){
			return {};
		}
		// return {};
		double d = c1.dist(c2);
		point p1 = {x1 + (x2 - x1) * r1 / d,
				   y1 + (y2 - y1) * r1 / d};
		point p2 = {x2 + (x1 - x2) * r2 / d,
					y2 + (y1 - y2) * r2 / d};

		p1 = c1.rotate(p1, PI/6.0);
		p2 = c2.rotate(p2, PI/6.0);
		
		line l1 = (line_segment){p1, p2}.get_line();
		line l2 = (line_segment){c1, c2}.get_line();

		vector<point> pts = l1.intersect(l2);
		vector<point> v1 = tangent_points(pts[0]);
		vector<point> v2 = c.tangent_points(pts[0]);

		if(l2.same_side(v1[0], v2[0])) swap(v1[0], v1[1]);

		vector<line_segment > ans;
		ans.push_back((line_segment){v1[0], v2[0]});
		ans.push_back((line_segment){v1[1], v2[1]});
		return ans;
	}

	// returns top tangents between the circles as vector of line segments
	// All cases handled. Cross tangent is not possible if one circle is completely inside the other.
	// In case of internal-touch, line segment of zero length passing through point of contact is returned.
	vector<line_segment> top_tangents(circle c){
		double x1 = xc, y1 = yc, r1 = r;
		double x2 = c.xc, y2 = c.yc, r2 = c.r;

		point c1 = {x1, y1};
		point c2 = {x2, y2};

		// one tangent : returns a zero length line segment through the point of contact 
		if(zero(c1.dist2(c2) - (r1-r2)*(r1-r2))){
			c1 = intersect(c)[0];
			return {{c1, c1}};
		}

		// circles inside
		if(c1.dist2(c2) + eps < (r1-r2)*(r1-r2)){
			return {};
		}

		if(abs(r1 - r2) <= eps){ // same radius
			line ll = (line_segment){c1, c2}.get_line();
			line l1 = ll.perp(c1);
			line l2 = ll.perp(c2);
			vector<point> v1 = l1.points_at_d(c1, r1);
			vector<point> v2 = l2.points_at_d(c2, r2);
			if(!zero(v1[0].dist2(v2[0]) - c1.dist2(c2))) {
				swap(v1[0], v1[1]);
			}
			vector<line_segment > ans;
			ans.push_back((line_segment){v1[0], v2[0]});
			ans.push_back((line_segment){v1[1], v2[1]});
			return ans;
		}

		double d = c1.dist(c2);
		point p1 = {x1 + (x2 - x1) * r1 / d,
				   y1 + (y2 - y1) * r1 / d};
		point p2 = {x2 + (x1 - x2) * r2 / d,
					y2 + (y1 - y2) * r2 / d};
		
		p1 = c1.rotate(p1, PI/6.0);
		p2 = c2.rotate(p2, 7.0*PI/6.0);

		line l1 = (line_segment){p1, p2}.get_line();
		line l2 = (line_segment){c1, c2}.get_line();

		vector<point> pts = l1.intersect(l2);
		vector<point> v1 = tangent_points(pts[0]);
		vector<point> v2 = c.tangent_points(pts[0]);

		if(!l2.same_side(v1[0], v2[0])) swap(v1[0], v1[1]);

		vector<line_segment > ans;
		ans.push_back((line_segment){v1[0], v2[0]});
		ans.push_back((line_segment){v1[1], v2[1]});
		return ans;
	}
};
