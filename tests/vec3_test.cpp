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
	CHECK(x==x+y-y);
	x+=y;
	x-=y;
	CHECK(x==cp);
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

TEST_CASE("propper arithmetic for vec3"){
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
	fp::vec3<double> zero_vec{0, 0, 0};
	CHECK(x.cross(x)==zero_vec);
	}
  SECTION("property test: antisymmetry of cross product"){
	constexpr int numTrials = 3;
	constexpr int min=-1.e5, max=1.e5;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> y{vec[3], vec[4], vec[5]};
	CHECK(x.cross(y)==(-1)*y.cross(x));
	}
  SECTION("property test: cross product is orthogonal to crossed vectors"){
	constexpr int numTrials = 3;
	constexpr double min=-100., max=100.;
	const std::vector<double> vec = GENERATE(take(numTrials,chunk<double>(6,random<double>(min,max))));
	fp::vec3<double> x{vec[0], vec[1], vec[2]};
	fp::vec3<double> y{vec[3], vec[4], vec[5]};


	fp::vec3<double> z=x.cross(y);
	auto target = Approx(0).margin(EPSILON);
	CHECK(z.dot(x)==target);
	CHECK(z.dot(y)==target);
  }

}

TEST_CASE("check -v correctness"){
    auto v = fp::vec3<float>{1.3, 6.8, 2.4};
    auto v_min = fp::vec3<float>{-1.3, -6.8, -2.4};
    CHECK(-v==v_min);

}
