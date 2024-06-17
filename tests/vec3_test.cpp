#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "flippy.hpp"


using Approx = Catch::Approx;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

float SMALL_FRACTION = 1e-5f;
const auto ZERO_ISH = WithinAbs(0, 0.00001);
auto ish(auto num){
    return (WithinRel(static_cast<float>(num),
                     static_cast<float> (SMALL_FRACTION)));
//    ||
//            WithinAbs(0, 0.000001));
}

auto operator"" _ish(long double value) { return ish(value); }
auto operator"" _ish(unsigned long long value) { return ish(value); }

TEMPLATE_TEST_CASE("proper initiation for vec3", "", float, double ){
    using REAL = TestType;


    SECTION("instantiation values are correct") {
        fp::vec3<REAL> v0{1, 12, 3};

        CHECK_THAT(v0[0], 1_ish);
        CHECK_THAT(v0[1], 12_ish);
        CHECK_THAT(v0[2], 3_ish);
    }

    SECTION("modify post instantiation"){
        fp::vec3<REAL> v0{1.2f, 4.f, 3.f};
        CHECK_THAT(v0[0], 1.2_ish);
        v0[0]+=1.1f;
        CHECK_THAT(v0[0], 2.3_ish);

    }

    fp::vec3<REAL> v0{1.2f, 4.f, 3.f};
    SECTION("make new vec from old"){
        fp::vec3<REAL> v1(v0);
        CHECK(v1==v0);
    }
    SECTION("make new vec from old"){
        fp::vec3<REAL> v1=v0;
        CHECK(v1==v0);
    }

    SECTION("check that copying works"){
        fp::vec3<REAL> v1{};
        v1=v0;
        CHECK(v1==v0);

    }

}

TEMPLATE_TEST_CASE("member function and associated operator checks", "",
                   float, double){
    using REAL = TestType;

    SECTION("case tests: add, 1"){
        fp::vec3<REAL> v0{1., 8., 17.};
        fp::vec3<REAL> v1{0, 1, 1};
        fp::vec3<REAL> sum{1., 9., 18.};
        CHECK(v0+v1==sum);
        v0+=v1;
        CHECK(v0==sum);
    }

    SECTION("case tests: add, 2"){
        fp::vec3<REAL> v0{1., 8., -17.};
        fp::vec3<REAL> v1{0, 1, 1};
        fp::vec3<REAL> sum{1., 9., -16.};
        CHECK(v0+v1==sum);
        v0+=v1;
        CHECK(v0==sum);
    }

    SECTION("case tests: subtract, 1"){
        fp::vec3<REAL> v0{12.1f, 3., -17.};
        fp::vec3<REAL> v1{0.2f, 2, 6};
        fp::vec3<REAL> diff{11.9f, 1., -23.};
        auto res = v0-v1-diff;
        CHECK_THAT(res.x, ZERO_ISH);
        CHECK_THAT(res.y, ZERO_ISH);
        CHECK_THAT(res.z, ZERO_ISH);
        v0-=v1;
        res = v0-diff;

        CHECK_THAT(res[0], ZERO_ISH);
        CHECK_THAT(res[1], ZERO_ISH);
        CHECK_THAT(res[2], ZERO_ISH);
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
        const std::vector<REAL> vec = GENERATE(take(numTrials,chunk(6,random<REAL>(min,max))));
        fp::vec3<REAL> x{vec[0], vec[1], vec[2]};
        fp::vec3<REAL> cp = x;
        fp::vec3<REAL> y{vec[3], vec[4], vec[5]};

        CHECK_THAT( (x+y-y).x, ish(x.x) );
        CHECK_THAT( (x+y-y).y, ish(x.y) );
        CHECK_THAT( (x+y-y).z, ish(x.z) );
        x+=y;
        x-=y;
        CHECK_THAT( cp.x, ish(x.x) );
        CHECK_THAT( cp.y, ish(x.y) );
        CHECK_THAT( cp.z, ish(x.z) );
    }

    SECTION("case test: norm"){
        fp::vec3<REAL> v{3., 4., 5.};
        CHECK_THAT(v.norm(), ish(sqrt(50)));
    }

    SECTION("case test: normalize"){
        fp::vec3<REAL> v{3., 4., 5.};
        v.normalize();
        CHECK_THAT(v.norm(), 1_ish);
    }

    SECTION("case test: normalize"){
        fp::vec3<double> v{3., 4., 5.};
        fp::vec3<double> v_norm = v/v.norm();
        v.normalize();
        CHECK_THAT(v.x, ish(v_norm.x));
        CHECK_THAT(v.y, ish(v_norm.y));
        CHECK_THAT(v.z, ish(v_norm.z));
    }


    SECTION("property test: norm square is self dot"){
        constexpr int numTrials = 3;
        constexpr int min=-1.e5, max=1.e5;
        const std::vector<REAL> vec = GENERATE(take(numTrials,chunk(3,random<REAL> (min,max))));
        fp::vec3<double> x{vec[0], vec[1], vec[2]};
        fp::floating_point_number auto x_norm = x.norm();
        CHECK_THAT(x.dot(x), ish(x_norm*x_norm));
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
        CHECK_THAT(cp.x, ZERO_ISH);
        CHECK_THAT(cp.y, ZERO_ISH);
        CHECK_THAT(cp.z, ZERO_ISH);
    }

    SECTION("property test: antisymmetry of cross product"){
        constexpr int numTrials = 3;
        constexpr int min=-1.e5, max=1.e5;
        const std::vector<double> vec = GENERATE(take(numTrials,chunk(6,random<double>(min,max))));
        fp::vec3<double> x{vec[0], vec[1], vec[2]};
        fp::vec3<double> y{vec[3], vec[4], vec[5]};
        fp::vec3<double> diff = x.cross(y) + y.cross(x);
        CHECK_THAT(diff.x, ZERO_ISH);
        CHECK_THAT(diff.y, ZERO_ISH);
        CHECK_THAT(diff.z, ZERO_ISH);
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
