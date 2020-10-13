#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen\Dense>

#pragma warning(disable:4018)
#pragma warning(disable:4244)


vec3 AInterpolatorVec3::LERP(vec3 a, vec3 b, double u)
{
	return a * (1 - u) + b * u;
}

// Reference: https://stackoverflow.com/questions/28036652/finding-the-shortest-distance-between-two-angles/28037434
vec3 AInterpolatorVec3::ShortestPathEuler(vec3 angle0, vec3 angle1)
{
	double xAngle = std::fmod(angle1[0] - angle0[0] + 180.0, 360.0) - 180.0;
	xAngle = xAngle < -180.0 ? xAngle + 360.0 : xAngle;

	double yAngle = std::fmod(angle1[1] - angle0[1] + 180.0, 360.0) - 180.0;
	yAngle = yAngle < -180.0 ? yAngle + 360.0 : yAngle;

	double zAngle = std::fmod(angle1[2] - angle0[2] + 180.0, 360.0) - 180.0;
	zAngle = zAngle < -180.0 ? zAngle + 360.0 : zAngle;

	return vec3(xAngle, yAngle, zAngle);
}

ASplineVec3::ASplineVec3() : mInterpolator(new ABernsteinInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
    if (mInterpolator) delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
    mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
    return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
    mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
    return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
    double fps = getFramerate();

	if (mInterpolator) { delete mInterpolator; }
    switch (type)
    {
	case LINEAR: mInterpolator = new ALinearInterpolatorVec3(); break;
	case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3(); break;
	case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3(); break;
	case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3(); break;
	case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3(); break;
	case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3(); break;
	case LINEAR_EULER: mInterpolator = new AEulerLinearInterpolatorVec3(); break;
	case CUBIC_EULER: mInterpolator = new AEulerCubicInterpolatorVec3(); break;
    };
    
    mInterpolator->setFramerate(fps);
    computeControlPoints();
    cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
    return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys[keyID].second = value;
    computeControlPoints();
    cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0)
    {
        mStartPoint = value;
        computeControlPoints(false);
    }
    else if (ID == mCtrlPoints.size() + 1)
    {
        mEndPoint = value;
		computeControlPoints(false);
    }
    else mCtrlPoints[ID-1] = value;
    cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
    mKeys.push_back(Key(time, value));

    if (updateCurve)
    {
        computeControlPoints();
        cacheCurve();
    }
}

int ASplineVec3::insertKey(double time, const vec3& value, bool updateCurve)
{
	if (mKeys.size() == 0)
	{
		appendKey(time, value, updateCurve);
		return 0;
	}

	for (int i = 0; i < mKeys.size(); ++i)
	{
		assert(time != mKeys[i].first);
		if (time < mKeys[i].first)
		{
			mKeys.insert(mKeys.begin() + i, Key(time, value));
			if (updateCurve)
			{
				computeControlPoints();
				cacheCurve();
			}
			return i;
		}
	}

	// Append at the end of the curve
	appendKey(time, value, updateCurve);
	return mKeys.size() - 1;
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
    if (mKeys.size() == 0)
    {
        appendKey(0, value, updateCurve);
    }
    else
    {
        double lastT = mKeys[mKeys.size() - 1].first;
        appendKey(lastT + 1, value, updateCurve);
    }
}

void ASplineVec3::deleteKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys.erase(mKeys.begin() + keyID);
    computeControlPoints();
    cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID) const
{
    assert(keyID >= 0 && keyID < mKeys.size());
    return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
    return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID) const
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0) return mStartPoint;
    else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
    else return mCtrlPoints[ID-1];
}

int ASplineVec3::getNumControlPoints() const
{
    return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
    mKeys.clear();
}

double ASplineVec3::getDuration() const 
{
    return mKeys.size() == 0 ? 0 : mKeys[mKeys.size()-1].first;
}

double ASplineVec3::getNormalizedTime(double t) const 
{
    return (t / getDuration());
}

double ASplineVec3::getKeyTime(int keyID) const
{
	assert(keyID >= 0 && keyID < mKeys.size());
	return mKeys[keyID].first;
}

vec3 ASplineVec3::getValue(double t) const
{
    if (mCachedCurve.size() == 0 || mKeys.size() == 0) return vec3();
	if (t < mKeys[0].first)
		return mCachedCurve[0];
	else
		t -= mKeys[0].first;

    double dt = mInterpolator->getDeltaTime();
    int rawi = (int)(t / dt); // assumes uniform spacing
    double frac = (t - rawi*dt) / dt;

	int i = mLooping? rawi % mCachedCurve.size() : std::min<int>(rawi, mCachedCurve.size() - 1);
	int inext = mLooping ? (i + 1) % mCachedCurve.size() : std::min<int>(i + 1, mCachedCurve.size() - 1);

    vec3 v1 = mCachedCurve[i];
    vec3 v2 = mCachedCurve[inext];
    vec3 v = v1*(1 - frac) + v2 * frac;
    return v;
}

void ASplineVec3::cacheCurve()
{
    mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
	if (mInterpolator->getType() == ASplineVec3::CUBIC_HERMITE)
	{
		std::vector<vec3> reworkCtrlPoints = std::vector<vec3>();
		for (int i = 2; i < mCtrlPoints.size(); i += 4)
		{
			reworkCtrlPoints.push_back(mCtrlPoints[i]);
		}
		reworkCtrlPoints.push_back(mCtrlPoints[mCtrlPoints.size() - 1]);
		mCtrlPoints.clear();
		for (int i = 0; i < reworkCtrlPoints.size(); ++i)
		{
			mCtrlPoints.push_back(reworkCtrlPoints[i]);
		}
	}
}

void ASplineVec3::computeControlPoints(bool updateEndPoints)
{
	if (mKeys.size() >= 2 && updateEndPoints)
	{
		int totalPoints = mKeys.size();

		//If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
		//They lie on the tangent of the first and last interpolation points.
		vec3 tmp = mKeys[0].second - mKeys[1].second;
		double n = tmp.Length();
		mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25; // distance to endpoint is 25% of distance between first 2 points

		tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
		n = tmp.Length();
		mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
	}
    mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

vec3* ASplineVec3::getCachedCurveData()
{
	return mCachedCurve.data();
}

vec3 * ASplineVec3::getControlPointsData()
{
	return mCtrlPoints.data();
}

int ASplineVec3::getNumCurveSegments() const
{
    return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
    return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
    mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
    return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
    return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys, 
    const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0;

	curve.clear();

	int numSegments = keys.size() - 1;
	for (int segment = 0; segment < numSegments; segment++)
    {
        for (double t = keys[segment].first; t < keys[segment+1].first - FLT_EPSILON; t += mDt)
        {
			// TODO: Compute u, fraction of duration between segment and segmentnext, for example,
			// u = 0.0 when t = keys[segment-1].first  
			// u = 1.0 when t = keys[segment].first
			u = t - keys[segment].first;

            val = interpolateSegment(keys, ctrlPoints, segment, u);
            curve.push_back(val);
        }
    }
	// add last point
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
		curve.push_back(val);
	}
}


// Interpolate p0 and p1 so that t = 0 returns p0 and t = 1 returns p1
vec3 ALinearInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
	curveValue = LERP(key0, key1, u);
	return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	curveValue = (b0 * (1 - u) * (1 - u) * (1 - u)) +
				 (b1 * 3.f * u * (1 - u) * (1 - u)) +
				 (b2 * 3.f * u * u * (1 - u)) +
				 (b3 * u * u * u);

	return curveValue;
}

vec3 ACasteljauInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	vec3 b01 = LERP(b0, b1, u);
	vec3 b12 = LERP(b1, b2, u);
	vec3 b23 = LERP(b2, b3, u);

	vec3 b01_12 = LERP(b01, b12, u);
	vec3 b12_23 = LERP(b12, b23, u);

	curveValue = LERP(b01_12, b12_23, u);

	return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
	// Hint: Using Eigen::MatrixXd data representations for a matrix operations

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	Eigen::MatrixXd G(3, 4);
	G << b0[0], b1[0], b2[0], b3[0],
		 b0[1], b1[1], b2[1], b3[1],
		 b0[2], b1[2], b2[2], b3[2];

	Eigen::MatrixXd M(4, 4);
	M << 1.0, -3.0, 3.0, -1.0,
		 0.0, 3.0, -6.0, 3.0,
		 0.0, 0.0, 3.0, -3.0,
		 0.0, 0.0, 0.0, 1.0;

	Eigen::MatrixXd U(4, 1);
	U << 1.0, u, u*u, u*u*u;

	Eigen::MatrixXd F = G * M * U;

	curveValue = vec3(F(0, 0), F(1, 0), F(2, 0));

	return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 p0;
	vec3 p1;
	vec3 q0; // slope at p0
	vec3 q1; // slope at p1
	vec3 curveValue(0, 0, 0);

	// TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  

	p0 = ctrlPoints[4 * segment];
	p1 = ctrlPoints[4 * segment + 1];
	q0 = ctrlPoints[4 * segment + 2];
	q1 = ctrlPoints[4 * segment + 3];

	curveValue = p0 * (2.0 * u * u * u - 3.0 * u * u + 1.0) +
				 p1 * (-2.0 * u * u * u + 3.0 * u * u) +
				 q0 * (u * u * u - 2.0 * u * u + u) +
				 q1 * (u * u * u - u * u);

	return curveValue;
}

vec3 ABSplineInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 curveValue(0, 0, 0);

	// Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = curve interval on knot vector in which to interpolate
	//     t = time value	

	// Step 1: determine the index j
	// Step 2: compute the n nonzero Bspline Basis functions N given j
	// Step 3: get the corresponding control points from the ctrlPoints vector
	// Step 4: compute the Bspline curveValue at time t

	return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys, 
    std::vector<vec3>& ctrlPoints, 
    vec3& startPoint, vec3& endPoint)
{
    ctrlPoints.clear();
    if (keys.size() <= 1) return;

    for (int i = 1; i < keys.size(); i++)
    {
        vec3 b0, b1, b2, b3;
		// TODO: compute b0, b1, b2, b3

		vec3 s0, s1;

		if (i != 1)
		{
			s0 = (keys[i].second - keys[i - 2].second) * 0.5f;
		}
		else
		{
			s0 = (keys[i].second - startPoint) * 0.5f;
		}
		if (i != keys.size() - 1)
		{
			s1 = (keys[i + 1].second - keys[i - 1].second) * 0.5f;
		}
		else
		{
			s1 = (endPoint - keys[i - 1].second) * 0.5f;
		}

		b0 = keys[i - 1].second;
		b3 = keys[i].second;

		b1 = b0 + s0 / 3;
		b2 = b3 - s1 / 3;

        ctrlPoints.push_back(b0);
        ctrlPoints.push_back(b1);
        ctrlPoints.push_back(b2);
        ctrlPoints.push_back(b3);
    }
}

void AHermiteInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys,
    std::vector<vec3>& ctrlPoints,
    vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	//ctrlPoints.resize(keys.size(), vec3(0, 0, 0));
	if (keys.size() <= 1) return;

	// TODO: 
	// For each key point pi, compute the corresonding value of the slope pi_prime.
	// Hints: Using Eigen::MatrixXd for a matrix data structures, 
	// this can be accomplished by solving the system of equations AC=D for C.
	// Don't forget to save the values computed for C in ctrlPoints
	// For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
	// For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

	// Step 1: Initialize A
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(keys.size(), keys.size());
	A(0, 0) = 2.0;
	A(0, 1) = 1.0;
	A(keys.size() - 1, keys.size() - 1) = 2.0;
	A(keys.size() - 1, keys.size() - 2) = 1.0;

	for (int i = 1; i < keys.size() - 1; ++i)
	{
		A(i, i - 1) = 1.0;
		A(i, i) = 4.0;
		A(i, i + 1) = 1.0;
	}

	//std::cout << "Here is the matrix A:\n" << A << std::endl;

	// Step 2: Initialize D
	Eigen::MatrixXd D = Eigen::MatrixXd(keys.size(), 3);
	D(0, 0) = 3.0 * (keys[1].second[0] - keys[0].second[0]);
	D(0, 1) = 3.0 * (keys[1].second[1] - keys[0].second[1]);
	D(0, 2) = 3.0 * (keys[1].second[2] - keys[0].second[2]);
	D(keys.size() - 1, 0) = 3.0 * (keys[keys.size() - 1].second[0] - keys[keys.size() - 2].second[0]);
	D(keys.size() - 1, 1) = 3.0 * (keys[keys.size() - 1].second[1] - keys[keys.size() - 2].second[1]);
	D(keys.size() - 1, 2) = 3.0 * (keys[keys.size() - 1].second[2] - keys[keys.size() - 2].second[2]);
	for (int i = 1; i < keys.size() - 1; ++i)
	{
		D(i, 0) = 3.0 * (keys[i + 1].second[0] - keys[i - 1].second[0]);
		D(i, 1) = 3.0 * (keys[i + 1].second[1] - keys[i - 1].second[1]);
		D(i, 2) = 3.0 * (keys[i + 1].second[2] - keys[i - 1].second[2]);
	}

	//std::cout << "Here is the matrix D:\n" << D << std::endl;

	// Step 3: Solve AC=D for C

	Eigen::MatrixXd C = A.inverse() * D;

	//std::cout << "Here is the matrix C:\n" << C << std::endl;

	// Step 4: Save control points in ctrlPoints

	for (int i = 0; i < keys.size() - 1; ++i)
	{
		ctrlPoints.push_back(keys[i].second);
		ctrlPoints.push_back(keys[i+1].second);
		ctrlPoints.push_back(vec3(C(i, 0), C(i, 1), C(i, 2)));
		ctrlPoints.push_back(vec3(C(i+1, 0), C(i+1, 1), C(i+1, 2)));
	}
}

void ABSplineInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys,
    std::vector<vec3>& ctrlPoints, 
    vec3& startPt, vec3& endPt)
{
    ctrlPoints.clear();
	ctrlPoints.resize(keys.size() + 2, vec3(0, 0, 0));
    if (keys.size() <= 1) return;

	// TODO:
	// Hints: 
	// 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

	// 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = interval on knot vector in which to interpolate
	//     t = time value
	//     l = derivative (l = 1 => 1st derivative)

	// Step 1: Calculate knot vector using a uniform BSpline
	//         (assune knots are evenly spaced 1 apart and the start knot is at time = 0.0)

	// Step 2: Calculate A matrix  for a natural BSpline
	//         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

	// Step 3: Calculate  D matrix composed of our target points to interpolate

	// Step 4: Solve AC=D for C 

	// Step 5: save control points in ctrlPoints
}


vec3 AEulerLinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys, 
	const std::vector<vec3>& ctrlPoints, 
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO:
	// Linear interpolate between key0 and key1
	// You should convert the angles to find the shortest path for interpolation

	key1 = key0 + ShortestPathEuler(key0, key1);

	curveValue = LERP(key0, key1, u);

	return curveValue;
}

vec3 AEulerCubicInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys, 
	const std::vector<vec3>& ctrlPoints, int segment, double t)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	// You should convert the angles to find the shortest path for interpolation

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	curveValue = (b0 * (1 - t) * (1 - t) * (1 - t)) +
		(b1 * 3.f * t * (1 - t) * (1 - t)) +
		(b2 * 3.f * t * t * (1 - t)) +
		(b3 * t * t * t);

	return curveValue;
}

void AEulerCubicInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys, 
	std::vector<vec3>& ctrlPoints, vec3 & startPoint, vec3 & endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	// Hint: One naive way is to first convert the keys such that the differences of the x, y, z Euler angles 
	//		 between every two adjacent keys are less than 180 degrees respectively 

	std::vector<vec3> convertedKeys = std::vector<vec3>();
	convertedKeys.push_back(keys[0].second);

	for (int i = 1; i < keys.size(); i++)
	{
		convertedKeys.push_back(convertedKeys[i - 1] + ShortestPathEuler(convertedKeys[i - 1], keys[i].second));
	}

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 b0, b1, b2, b3;

		// TODO: compute b0, b1, b2, b3

		vec3 s0, s1;

		if (i != 1)
		{
			s0 = (convertedKeys[i] - convertedKeys[i - 2]) * 0.5f;
		}
		else
		{
			s0 = convertedKeys[i] - convertedKeys[i - 1];
		}
		if (i != keys.size() - 1)
		{
			s1 = (convertedKeys[i + 1] - convertedKeys[i - 1]) * 0.5f;
		}
		else
		{
			s1 = convertedKeys[i] - convertedKeys[i - 1];
		}

		b0 = convertedKeys[i - 1];
		b3 = convertedKeys[i];

		b1 = b0 + s0 / 3;
		b2 = b3 - s1 / 3;

		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
}
