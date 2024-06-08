#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "flippy.hpp"

const double EPSILON = 1e-9;

TEST_CASE("proper initiation for vec3"){

  SECTION("instantiation values are correct") {
	fp::vec3<double> v0{1, 12, 3};
	REQUIRE(v0[0]==1);
	REQUIRE(v0[1]==12);
	REQUIRE(v0[2]==3);
  }

  SECTION("modify post instantiation"){
	fp::vec3<double> v0{1.2, 4, 3};
	REQUIRE(v0[0]==1.2);
	v0[0]+=1.1;
	REQUIRE(v0[0]==2.3);

  }

  SECTION("make new vec from old"){
	fp::vec3<double> v0{1.2, 4, 3};
	fp::vec3<double> v1(v0);
	CHECK(v1==v0);
  }
  SECTION("make new vec from old"){
	fp::vec3<double> v0{1.2, 4, 3};
	fp::vec3<double> v1=v0;
	CHECK(v1==v0);
  }

  SECTION("check that copying works"){
	fp::vec3<double> v0{1.2, 4, 3};
	fp::vec3<double> v1{};
	v1=v0;
	CHECK(v1==v0);

  }

}

TEST_CASE("member function and associated operator checks"){
  SECTION("case tests: add, 1"){
    fp::vec3<double> v0{1., 8., 17.};
	fp::vec3<double> v1{0, 1, 1};
	fp::vec3<double> sum{1., 9., 18.};
	CHECK(v0+v1==sum);
	v0+=v1;
	CHECK(v0==sum);
  }

  SECTION("case tests: add, 2"){
	fp::vec3<double> v0{1., 8., -17.};
	fp::vec3<double> v1{0, 1, 1};
	fp::vec3<double> sum{1., 9., -16.};
	CHECK(v0+v1==sum);
	v0+=v1;
	CHECK(v0==sum);
  }

  SECTION("case tests: subtract, 1"){
	fp::vec3<double> v0{12.1, 3., -17.};
	fp::vec3<double> v1{0.2, 2, 6};
	fp::vec3<double> sum{11.9, 1., -23.};
	CHECK(v0-v1==sum);
	v0-=v1;
	CHECK(v0==sum);
  }

  SECTION("property test: in place addition/subtraction is the same as adding to itself"){
	constexpr int numTrials = 3;
	constexpr int min=-1.e5, max=1.e5;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> y{vec[3], vec[4], vec[5]};
	fp::vec3<double> cp = x;
	cp = cp + y;
	x += y;
	CHECK(x==cp);
	cp -= y;
	x = x - y;
	CHECK(x==cp);
  }

  SECTION("property test: ad and subtract cancel"){
	constexpr int numTrials = 3;
	constexpr int min=-1.e5, max=1.e5;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> cp = x;
	fp::vec3<double> y{vec[3], vec[4], vec[5]};
    auto target = Approx(0).margin(EPSILON);
	CHECK((x-(x+y-y)).x==target);
	CHECK((x-(x+y-y)).y==target);
	CHECK((x-(x+y-y)).z==target);
	x+=y;
	x-=y;
	CHECK((x-cp).x==target);
	CHECK((x-cp).y==target);
	CHECK((x-cp).z==target);
  }

  SECTION("case test: norm"){
    fp::vec3<double> v{3., 4., 5.};
    auto target = Approx(sqrt(50)).margin(EPSILON);
    CHECK(v.norm()==target);
  }

    SECTION("case test: normalize"){
        fp::vec3<double> v{3., 4., 5.};
        auto target = Approx(sqrt(1)).margin(EPSILON);
        v.normalize();
        CHECK(v.norm()==target);
    }

    SECTION("case test: normalize"){
        fp::vec3<double> v{3., 4., 5.};
        fp::vec3<double> v_norm = v/v.norm();
        auto target = Approx(sqrt(0)).margin(EPSILON);
        v.normalize();
        CHECK((v.x-v_norm.x)==target);
        CHECK((v.y-v_norm.y)==target);
        CHECK((v.z-v_norm.z)==target);
    }


  SECTION("property test: norm square is self dot"){
	constexpr int numTrials = 3;
	constexpr int min=-1.e5, max=1.e5;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk(3,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	auto x_norm = x.norm();
	auto x_norm_square_ish = Approx(x_norm*x_norm).margin(EPSILON);
	CHECK(x.dot(x)==x_norm_square_ish);
  }

}

TEST_CASE("Operator checks"){
  SECTION("case test: Check that == operator works"){
	fp::vec3<double> x{1, 0, 0};
	fp::vec3<double> x_other{1, 0, 0};
	fp::vec3<double> y{0, 1, 0};
	CHECK_FALSE(x==y);
	CHECK(x==x);
	CHECK(x==x_other);
  }

}

TEST_CASE("proper arithmetic for vec3"){
    auto zero_ish = Approx(0).margin(1e-6);

	SECTION("test cross product x cross y is z"){
	fp::vec3<double> x{1, 0, 0};
	fp::vec3<double> y{0, 1, 0};
	fp::vec3<double> z{0, 0, 1};
	auto cp = x.cross(y);
	CHECK(z==cp);
  }
  SECTION("property test: self cross is zero"){
	constexpr int numTrials = 10;
	constexpr double min=-1e5, max=1e5;
	auto vec = GENERATE(take(numTrials,chunk(3,random(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
    auto cp = x.cross(x);
    CHECK(cp.x==zero_ish);
    CHECK(cp.y==zero_ish);
    CHECK(cp.z==zero_ish);
	}

  SECTION("property test: antisymmetry of cross product"){
	constexpr int numTrials = 3;
	constexpr int min=-1.e5, max=1.e5;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> y{vec[3], vec[4], vec[5]};
	fp::vec3<double> diff = x.cross(y) + y.cross(x);
    CHECK(diff.x==zero_ish);
    CHECK(diff.y==zero_ish);
    CHECK(diff.z==zero_ish);
//    CHECK(x.cross(y)==(-1)*y.cross(x));
	}
  SECTION("property test: cross product is orthogonal to crossed vectors"){
	constexpr int numTrials = 3;
	constexpr double min=-100., max=100.;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk<double>(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> y{vec[3], vec[4], vec[5]};


	fp::vec3<double> z=x.cross(y);
	CHECK(z.dot(x)==zero_ish);
	CHECK(z.dot(y)==zero_ish);
  }

}

TEST_CASE("check -v correctness"){
    SECTION("arithmetic correctness") {
        auto v = fp::vec3<float>{1.3f, 6.8f, 2.4f};
        auto v_min = fp::vec3<float>{-1.3f, -6.8f, -2.4f};
        CHECK(-v == v_min);
    }
    SECTION("rvalue - is correctly returned"){
//        auto make_min = [](fp::vec3<float>&& v){return -v;};
//        auto v_min = fp::vec3<float>{-1.3f, -6.8f, -2.4f};
        auto v = fp::vec3<float>{1.3f, 6.8f, 2.4f};
        CHECK(v ==  -(fp::vec3<float>{-1.3f, -6.8f, -2.4f}));
    }

    SECTION("rvalue - is correctly returned 2"){
        auto make_min = [](fp::vec3<float>&& v){return -v;};
        auto v_min = fp::vec3<float>{-1.3f, -6.8f, -2.4f};
        fp::vec3<float> temp{1.3f, 6.8f, 2.4f};
        auto v = make_min(std::move(temp));
        CHECK(v == v_min);
    }

}
