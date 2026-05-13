
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using std::cout;
using std::fixed;
using std::max;
using std::min;
using std::numeric_limits;
using std::setprecision;
using std::size_t;
using std::string;
using std::vector;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kLakeRadius = 10.0;
constexpr double kHareSpeed = 1.0;
constexpr double kCoyoteSpeed = 3.5;
constexpr double kRequiredLeadTime = 0.2;
constexpr int kMaxDirectionChanges = 5;
constexpr double kMaxAngleGap = 179.0;
constexpr double kEps = 1e-12;

struct Point2D {
    double x{0.0};
    double y{0.0};
};

struct Node {
    int id{-1};
    Point2D position{};
};

struct Edge {
    int from{-1};
    int to{-1};
};

struct ValidationResult {
    bool ok{false};
    string reason;
    double hareArrivalTime{0.0};
    double winningMargin{0.0};
    Point2D shorePoint{};
    int directionChanges{0};
};

class Geometry {
public:
    static double distance(const Point2D& a, const Point2D& b) {
        return std::hypot(a.x - b.x, a.y - b.y);
    }

    static bool nearlyZero(double value) {
        return std::fabs(value) < kEps;
    }

    static double normalizeDegrees(double angle) {
        angle = std::fmod(angle, 360.0);
        if (angle < 0.0) angle += 360.0;
        return angle;
    }

    static double pointAngleDegrees(const Point2D& p) {
        if (nearlyZero(p.x) && nearlyZero(p.y)) {
            return numeric_limits<double>::quiet_NaN();
        }
        return normalizeDegrees(std::atan2(p.y, p.x) * 180.0 / kPi);
    }

    static bool insideLake(const Point2D& p) {
        return (p.x * p.x + p.y * p.y) <= (kLakeRadius * kLakeRadius + kEps);
    }

    static int directionId(const Point2D& a, const Point2D& b) {
        const int dx = static_cast<int>(std::llround(b.x - a.x));
        const int dy = static_cast<int>(std::llround(b.y - a.y));

        if (dx == 1 && dy == 0)   return 0;
        if (dx == 1 && dy == 1)   return 1;
        if (dx == 0 && dy == 1)   return 2;
        if (dx == -1 && dy == 1)  return 3;
        if (dx == -1 && dy == 0)  return 4;
        if (dx == -1 && dy == -1) return 5;
        if (dx == 0 && dy == -1)  return 6;
        if (dx == 1 && dy == -1)  return 7;
        return 8;
    }

    static int countDirectionChanges(const vector<Point2D>& path) {
        int previous = 8;
        int changes = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const int current = directionId(path[i - 1], path[i]);
            if (previous != 8 && current != previous) {
                ++changes;
            }
            previous = current;
        }
        return changes;
    }

    static bool intersectCircle(const Point2D& a, const Point2D& b, Point2D& intersectionPoint) {
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double A = dx * dx + dy * dy;
        const double B = 2.0 * (a.x * dx + a.y * dy);
        const double C = a.x * a.x + a.y * a.y - kLakeRadius * kLakeRadius;

        double discriminant = B * B - 4.0 * A * C;
        if (discriminant < -kEps) {
            return false;
        }

        discriminant = std::max(0.0, discriminant);
        const double sqrtD = std::sqrt(discriminant);
        const double t1 = (-B - sqrtD) / (2.0 * A);
        const double t2 = (-B + sqrtD) / (2.0 * A);

        double best = 2.0;
        bool found = false;
        for (double t : {t1, t2}) {
            if (t >= -kEps && t <= 1.0 + kEps && t < best) {
                best = t;
                found = true;
            }
        }

        if (!found) {
            return false;
        }

        const double clamped = std::min(std::max(best, 0.0), 1.0);
        intersectionPoint = Point2D{a.x + clamped * dx, a.y + clamped * dy};
        return true;
    }

    static vector<Point2D> expandPolyline(const vector<Point2D>& anchors) {
        vector<Point2D> expanded;
        if (anchors.empty()) return expanded;

        expanded.push_back(anchors.front());
        for (size_t i = 1; i < anchors.size(); ++i) {
            const Point2D a = anchors[i - 1];
            const Point2D b = anchors[i];
            const int dx = static_cast<int>(std::llround(b.x - a.x));
            const int dy = static_cast<int>(std::llround(b.y - a.y));
            const int steps = std::max(1, std::gcd(std::abs(dx), std::abs(dy)));

            const double stepX = static_cast<double>(dx) / steps;
            const double stepY = static_cast<double>(dy) / steps;

            for (int s = 1; s <= steps; ++s) {
                expanded.push_back(Point2D{a.x + stepX * s, a.y + stepY * s});
            }
        }
        return expanded;
    }

    static vector<Point2D> interiorRouteAnchors() {
        return {
            {0.0, 0.0},
            {0.0, 1.0},
            {-1.0, 1.0},
            {-1.0, -1.0},
            {0.0, -1.0}
        };
    }

    static Point2D outsideProbe() {
        return {7.0, -8.0};
    }
};

class HareCoyoteModel {
public:
    double lakeRadius() const { return kLakeRadius; }
    double hareSpeed() const { return kHareSpeed; }
    double coyoteSpeed() const { return kCoyoteSpeed; }
    double requiredLeadTime() const { return kRequiredLeadTime; }
    int maxDirectionChanges() const { return kMaxDirectionChanges; }

    double coyoteAngularSpeedDegPerSec() const {
        return (kCoyoteSpeed / kLakeRadius) * 180.0 / kPi;
    }

    double coyoteAngleDegrees(double elapsedTime) const {
        return Geometry::normalizeDegrees(coyoteAngularSpeedDegPerSec() * elapsedTime);
    }

    static bool angleConstraintSatisfied(double hareAngle, double coyoteAngle) {
        if (std::isnan(hareAngle)) return true;
        const double gap = Geometry::normalizeDegrees(hareAngle - coyoteAngle);
        return gap > 0.0 && gap < kMaxAngleGap;
    }

    static double angleGapDegrees(double hareAngle, double coyoteAngle) {
        return Geometry::normalizeDegrees(hareAngle - coyoteAngle);
    }

    double winningMarginSeconds(double hareAngle, double coyoteAngle) const {
        return angleGapDegrees(hareAngle, coyoteAngle) / coyoteAngularSpeedDegPerSec();
    }

    double coyoteArrivalTimeAtSameAngle(double hareArrivalTime, double winningMargin) const {
        return hareArrivalTime + winningMargin;
    }

    bool validatePolyline(const vector<Point2D>& detailedPath, double& elapsedTime) const {
        elapsedTime = 0.0;
        if (detailedPath.size() < 2) return false;

        for (size_t i = 1; i < detailedPath.size(); ++i) {
            const Point2D a = detailedPath[i - 1];
            const Point2D b = detailedPath[i];
            const double segmentLength = Geometry::distance(a, b);
            const int samples = max(25, min(120, static_cast<int>(std::ceil(segmentLength / 0.01))));

            for (int s = 0; s <= samples; ++s) {
                const double u = static_cast<double>(s) / samples;
                const Point2D current{
                    a.x + (b.x - a.x) * u,
                    a.y + (b.y - a.y) * u
                };

                if (!Geometry::insideLake(current) &&
                    !(Geometry::nearlyZero(current.x - b.x) && Geometry::nearlyZero(current.y - b.y))) {
                    return false;
                }

                if (Geometry::nearlyZero(current.x) && Geometry::nearlyZero(current.y)) {
                    continue;
                }

                const double currentTime = elapsedTime + (segmentLength * u / kHareSpeed);
                const double hareAngle = Geometry::pointAngleDegrees(current);
                const double coyoteAngle = coyoteAngleDegrees(currentTime);

                if (!angleConstraintSatisfied(hareAngle, coyoteAngle)) {
                    return false;
                }
            }

            elapsedTime += segmentLength / kHareSpeed;
        }

        return true;
    }

private:
};

class RouteGraph {
public:
    int addNode(const Point2D& position) {
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{id, position});
        return id;
    }

    void addEdge(int from, int to) {
        edges_.push_back(Edge{from, to});
    }

private:
    vector<Node> nodes_;
    vector<Edge> edges_;
};

class ProgramController {
public:
    void run() {
        const vector<Point2D> interiorAnchors = Geometry::interiorRouteAnchors();
        const Point2D probe = Geometry::outsideProbe();

        Point2D shorePoint{};
        if (!Geometry::intersectCircle(interiorAnchors.back(), probe, shorePoint)) {
            cout << "Validation failed for the built-in path.\n";
            cout << "Could not compute the final shore crossing.\n";
            return;
        }

        vector<Point2D> fullRouteAnchors = interiorAnchors;
        fullRouteAnchors.push_back(shorePoint);

        RouteGraph graph;
        vector<int> nodeIds;
        nodeIds.reserve(fullRouteAnchors.size());
        for (const Point2D& p : fullRouteAnchors) {
            nodeIds.push_back(graph.addNode(p));
        }
        for (size_t i = 1; i < nodeIds.size(); ++i) {
            graph.addEdge(nodeIds[i - 1], nodeIds[i]);
        }

        ValidationResult result = validateRoute(fullRouteAnchors, shorePoint);
        if (!result.ok) {
            cout << "Validation failed for the built-in path.\n";
            if (!result.reason.empty()) {
                cout << result.reason << '\n';
            }
            return;
        }

        printReport(result, fullRouteAnchors);
    }

private:
    ValidationResult validateRoute(const vector<Point2D>& routeAnchors, const Point2D& shorePoint) const {
        ValidationResult result;

        double elapsedTime = 0.0;
        if (!model_.validatePolyline(routeAnchors, elapsedTime)) {
            result.reason = "The route violates the angle or boundary constraints.";
            return result;
        }

        const double hareArrivalTime = elapsedTime;
        const double hareAngle = Geometry::pointAngleDegrees(shorePoint);
        const double coyoteAngle = model_.coyoteAngleDegrees(hareArrivalTime);
        const double winningMargin = model_.winningMarginSeconds(hareAngle, coyoteAngle);
        const int directionChanges = Geometry::countDirectionChanges(routeAnchors);

        if (!HareCoyoteModel::angleConstraintSatisfied(hareAngle, coyoteAngle)) {
            result.reason = "The final shore crossing does not preserve the required angle ordering.";
            return result;
        }
        if (winningMargin < model_.requiredLeadTime()) {
            result.reason = "The hare does not beat the coyote by at least 0.2 seconds.";
            return result;
        }
        if (directionChanges > model_.maxDirectionChanges()) {
            result.reason = "The route uses more than 5 direction changes.";
            return result;
        }

        result.ok = true;
        result.hareArrivalTime = hareArrivalTime;
        result.winningMargin = winningMargin;
        result.shorePoint = shorePoint;
        result.directionChanges = directionChanges;
        return result;
    }

    void printReport(
        const ValidationResult& result,
        const vector<Point2D>& fullRouteAnchors
    ) const {
        cout << fixed << setprecision(6);
        cout << "Valid escape path found.\n\n";
        cout << "Lake radius: " << model_.lakeRadius() << " m\n";
        cout << "Hare speed: " << model_.hareSpeed() << " m/s\n";
        cout << "Coyote speed: " << model_.coyoteSpeed() << " m/s\n";
        cout << "Winning margin: " << result.winningMargin << " s\n";
        cout << "Hare arrival time: " << result.hareArrivalTime << " s\n";
        cout << "Coyote arrival time at same angle: "
             << model_.coyoteArrivalTimeAtSameAngle(result.hareArrivalTime, result.winningMargin) << " s\n";
        cout << "Hare angle at shore: " << Geometry::pointAngleDegrees(result.shorePoint) << " deg\n";
        cout << "Coyote angle at shore: " << model_.coyoteAngleDegrees(result.hareArrivalTime) << " deg\n";
        cout << "Direction changes used: " << result.directionChanges << "\n";
        cout << "Shore boundary check: x^2 + y^2 = "
             << (result.shorePoint.x * result.shorePoint.x + result.shorePoint.y * result.shorePoint.y)
             << "\n\n";

        cout << "Route waypoints (5 straight segments):\n";
        for (size_t i = 0; i < fullRouteAnchors.size(); ++i) {
            cout << i << ": (" << fullRouteAnchors[i].x << ", " << fullRouteAnchors[i].y << ")\n";
        }

        cout << "\nNote: intermediate lattice points that lie on a straight segment are not extra waypoints.\n";
        cout << "Final shore crossing: (" << result.shorePoint.x << ", " << result.shorePoint.y << ")\n";
    }

    HareCoyoteModel model_;
};

} // namespace

int main() {
    ProgramController controller;
    controller.run();
    return 0;
}