function C_xi = C_xi_h4_n4(dxi1,dxi2,dxi3,dxi4,dxi5,dxi6,dxi7,dxi8,dxi9,dxi10,dxi11,dxi12,dxi13,dxi14,...
    mu1,mu2,mu3,mu4,...
    xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11,xi12,xi13,xi14)
%C_xi_h4_n4
%    C_xi = C_xi_h4_n4(DXI1,DXI2,DXI3,DXI4,DXI5,DXI6,DXI7,DXI8,DXI9,DXI10,DXI11,DXI12,DXI13,DXI14,MU1,MU2,MU3,MU4,XI2,XI3,XI4,XI5,XI6,XI7,XI8,XI9,XI10,XI11,XI12,XI13,XI14)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    27-Nov-2022 19:36:18

t2 = mu3+mu4;
t3 = dxi1.*xi10;
t4 = dxi1.*xi11;
t5 = dxi4.*xi10;
t6 = dxi1.*xi14;
t7 = dxi4.*xi11;
t8 = dxi5.*xi10;
t9 = dxi5.*xi11;
t10 = dxi4.*xi14;
t11 = dxi8.*xi10;
t12 = dxi5.*xi14;
t13 = dxi8.*xi11;
t14 = dxi9.*xi10;
t15 = dxi9.*xi11;
t16 = dxi8.*xi14;
t17 = dxi9.*xi14;
t18 = dxi12.*xi14;
t19 = dxi13.*xi14;
t20 = dxi14.*xi14;
t21 = mu2.*xi2;
t22 = mu2.*xi3;
t23 = mu3.*xi2;
t24 = mu3.*xi3;
t25 = mu4.*xi2;
t26 = mu4.*xi3;
t27 = mu2.*xi6;
t28 = mu3.*xi6;
t29 = mu3.*xi7;
t30 = mu4.*xi6;
t31 = mu4.*xi7;
t32 = mu3.*xi10;
t33 = mu4.*xi10;
t34 = mu4.*xi11;
t35 = xi4+xi5;
t36 = xi8+xi9;
t37 = xi12+xi13;
t41 = mu4.*xi14.*2.0;
t42 = dxi1+dxi4+dxi5;
t38 = t32.*2.0;
t39 = t33.*2.0;
t40 = t34.*2.0;
t43 = mu2+t2;
t44 = dxi1.*t21;
t45 = dxi1.*t22;
t46 = dxi1.*t23;
t47 = dxi1.*t24;
t48 = dxi1.*t25;
t49 = dxi1.*t26;
t50 = dxi1.*t27;
t51 = dxi1.*t28;
t52 = dxi1.*t29;
t53 = dxi1.*t30;
t54 = dxi1.*t31;
t55 = dxi4.*t27;
t56 = dxi4.*t28;
t57 = dxi5.*t27;
t58 = mu3.*t3;
t59 = dxi4.*t29;
t60 = dxi4.*t30;
t61 = dxi5.*t28;
t62 = dxi6.*t27;
t63 = mu4.*t3;
t64 = dxi4.*t31;
t65 = dxi5.*t29;
t66 = dxi5.*t30;
t67 = dxi6.*t28;
t68 = mu4.*t4;
t69 = dxi5.*t31;
t70 = dxi6.*t29;
t71 = dxi6.*t30;
t72 = dxi7.*t28;
t73 = mu3.*t5;
t74 = dxi6.*t31;
t75 = dxi7.*t29;
t76 = dxi7.*t30;
t77 = mu4.*t5;
t78 = mu3.*t8;
t79 = dxi7.*t31;
t80 = mu4.*t7;
t81 = mu4.*t8;
t82 = mu4.*t9;
t83 = mu3.*t11;
t84 = mu4.*t11;
t85 = mu3.*t14;
t86 = mu4.*t13;
t87 = mu4.*t14;
t88 = dxi10.*t32;
t89 = mu4.*t15;
t90 = dxi10.*t33;
t91 = dxi10.*t34;
t92 = dxi11.*t33;
t93 = dxi11.*t34;
t94 = mu4.*t20;
t95 = cos(t35);
t96 = cos(t36);
t97 = cos(t37);
t98 = sin(t35);
t99 = sin(t36);
t100 = sin(t37);
t101 = -t20;
t102 = t35+t36;
t103 = t36+t37;
t104 = dxi8+dxi9+t42;
t105 = dxi10.*t97;
t106 = dxi11.*t97;
t107 = dxi6.*t99;
t108 = dxi7.*t99;
t109 = dxi10.*t100;
t110 = dxi11.*t100;
t111 = dxi14.*t100;
t112 = t96.*xi6;
t113 = t96.*xi7;
t114 = t97.*xi10;
t115 = t97.*xi11;
t116 = t97.*xi14;
t117 = cos(t102);
t118 = cos(t103);
t119 = -t44;
t120 = -t45;
t121 = -t46;
t122 = -t47;
t123 = -t48;
t124 = -t49;
t125 = -t50;
t126 = -t51;
t127 = -t52;
t128 = -t53;
t129 = -t54;
t130 = -t55;
t131 = -t56;
t132 = -t57;
t133 = -t58;
t134 = -t59;
t135 = -t60;
t136 = -t61;
t137 = -t63;
t138 = -t64;
t139 = -t65;
t140 = -t66;
t141 = -t68;
t142 = -t69;
t143 = -t73;
t144 = -t77;
t145 = -t78;
t146 = -t80;
t147 = -t81;
t148 = -t82;
t149 = -t83;
t150 = -t84;
t151 = -t85;
t152 = -t86;
t153 = -t87;
t154 = -t89;
t155 = sin(t102);
t156 = sin(t103);
t157 = t37+t102;
t164 = t3.*t97;
t165 = t4.*t97;
t166 = t5.*t97;
t167 = t6.*t97;
t168 = t7.*t97;
t169 = t8.*t97;
t170 = t9.*t97;
t171 = t10.*t97;
t172 = t11.*t97;
t173 = t12.*t97;
t174 = t13.*t97;
t175 = t14.*t97;
t176 = t15.*t97;
t177 = t16.*t97;
t178 = t17.*t97;
t181 = t18.*t97;
t182 = t19.*t97;
t183 = dxi2.*mu3.*t98;
t184 = dxi2.*mu4.*t98;
t185 = dxi3.*mu3.*t98;
t186 = dxi3.*mu4.*t98;
t187 = dxi6.*mu2.*t98;
t188 = dxi6.*mu3.*t98;
t189 = dxi6.*mu4.*t98;
t190 = dxi7.*mu3.*t98;
t191 = dxi7.*mu4.*t98;
t196 = dxi10.*mu3.*t99;
t197 = dxi10.*mu4.*t99;
t198 = dxi11.*mu4.*t99;
t200 = t23.*t95;
t201 = t24.*t95;
t202 = t25.*t95;
t203 = t26.*t95;
t204 = t27.*t95;
t205 = t28.*t95;
t206 = t29.*t95;
t207 = t30.*t95;
t208 = t31.*t95;
t209 = t28.*t96;
t210 = t29.*t96;
t211 = t30.*t96;
t212 = t31.*t96;
t213 = t32.*t96;
t214 = t33.*t96;
t215 = t34.*t96;
t217 = t3.*t100;
t218 = t4.*t100;
t219 = t5.*t100;
t220 = t7.*t100;
t221 = t8.*t100;
t222 = t9.*t100;
t223 = t11.*t100;
t224 = t13.*t100;
t225 = t14.*t100;
t226 = t15.*t100;
t227 = dxi12+dxi13+t104;
t228 = t46.*t95;
t229 = t47.*t95;
t230 = t48.*t95;
t231 = t49.*t95;
t232 = t50.*t95;
t233 = t51.*t95;
t235 = dxi6.*t21.*t95;
t236 = t52.*t95;
t237 = t53.*t95;
t240 = dxi6.*t22.*t95;
t242 = t54.*t95;
t246 = t55.*t95;
t253 = t56.*t95;
t254 = t57.*t95;
t259 = t59.*t95;
t260 = t60.*t95;
t261 = t61.*t95;
t263 = t64.*t95;
t264 = t65.*t95;
t265 = t66.*t95;
t266 = t69.*t95;
t267 = t51.*t96;
t268 = t52.*t96;
t269 = t53.*t96;
t270 = t54.*t96;
t271 = t56.*t96;
t272 = t58.*t96;
t273 = t59.*t96;
t274 = t60.*t96;
t275 = t61.*t96;
t276 = t63.*t96;
t277 = t64.*t96;
t278 = t65.*t96;
t279 = t66.*t96;
t280 = t68.*t96;
t281 = t69.*t96;
t282 = t73.*t96;
t283 = t77.*t96;
t284 = t78.*t96;
t285 = t80.*t96;
t286 = t81.*t96;
t289 = t82.*t96;
t296 = t83.*t96;
t300 = t84.*t96;
t301 = t85.*t96;
t303 = t86.*t96;
t304 = t87.*t96;
t305 = t89.*t96;
t312 = dxi14.*t33.*t97;
t314 = dxi14.*t34.*t97;
t331 = t18.*t100.*xi10;
t332 = t18.*t100.*xi11;
t333 = t19.*t100.*xi10;
t334 = t19.*t100.*xi11;
t343 = t39.*t97;
t344 = t40.*t97;
t345 = t41.*t97;
t346 = t39.*t100.*xi14;
t347 = t40.*t100.*xi14;
t348 = dxi1.*t2.*t98;
t356 = t44.*t98.*xi6;
t357 = t45.*t98.*xi6;
t358 = t46.*t98.*xi6;
t359 = t46.*t98.*xi7;
t360 = t47.*t98.*xi6;
t361 = t48.*t98.*xi6;
t362 = t47.*t98.*xi7;
t363 = t48.*t98.*xi7;
t364 = t49.*t98.*xi6;
t365 = dxi4.*t21.*t98.*xi6;
t366 = t49.*t98.*xi7;
t367 = dxi4.*t22.*t98.*xi6;
t368 = dxi4.*t23.*t98.*xi6;
t369 = dxi5.*t21.*t98.*xi6;
t370 = dxi4.*t23.*t98.*xi7;
t371 = dxi4.*t24.*t98.*xi6;
t372 = dxi4.*t25.*t98.*xi6;
t373 = dxi5.*t22.*t98.*xi6;
t374 = dxi5.*t23.*t98.*xi6;
t375 = dxi4.*t24.*t98.*xi7;
t376 = dxi4.*t25.*t98.*xi7;
t377 = dxi4.*t26.*t98.*xi6;
t378 = dxi5.*t23.*t98.*xi7;
t379 = dxi5.*t24.*t98.*xi6;
t380 = dxi5.*t25.*t98.*xi6;
t381 = dxi4.*t26.*t98.*xi7;
t382 = dxi5.*t24.*t98.*xi7;
t383 = dxi5.*t25.*t98.*xi7;
t384 = dxi5.*t26.*t98.*xi6;
t385 = dxi5.*t26.*t98.*xi7;
t387 = t3.*t28.*t99;
t388 = t3.*t29.*t99;
t389 = t3.*t30.*t99;
t391 = t4.*t30.*t99;
t392 = t3.*t31.*t99;
t394 = t4.*t31.*t99;
t395 = t5.*t28.*t99;
t397 = t5.*t29.*t99;
t398 = t5.*t30.*t99;
t399 = t8.*t28.*t99;
t400 = t7.*t30.*t99;
t401 = t5.*t31.*t99;
t402 = t8.*t29.*t99;
t403 = t8.*t30.*t99;
t404 = t7.*t31.*t99;
t405 = t9.*t30.*t99;
t406 = t8.*t31.*t99;
t407 = t9.*t31.*t99;
t408 = t11.*t28.*t99;
t409 = t11.*t29.*t99;
t410 = t11.*t30.*t99;
t411 = t14.*t28.*t99;
t412 = t13.*t30.*t99;
t413 = t11.*t31.*t99;
t414 = t14.*t29.*t99;
t415 = t14.*t30.*t99;
t416 = t13.*t31.*t99;
t417 = t15.*t30.*t99;
t418 = t14.*t31.*t99;
t419 = t15.*t31.*t99;
t420 = t63.*t100.*xi14;
t421 = t68.*t100.*xi14;
t422 = t77.*t100.*xi14;
t423 = t80.*t100.*xi14;
t424 = t81.*t100.*xi14;
t425 = t82.*t100.*xi14;
t427 = t84.*t100.*xi14;
t428 = t86.*t100.*xi14;
t429 = t87.*t100.*xi14;
t430 = t89.*t100.*xi14;
t431 = t18.*t33.*t100;
t432 = t18.*t34.*t100;
t433 = t19.*t33.*t100;
t434 = t19.*t34.*t100;
t456 = mu4.*t42.*t99;
t457 = dxi1.*t43.*t98;
t716 = mu4.*t99.*t104;
t717 = mu4.*t100.*t104;
t718 = t2.*t42.*t98;
t719 = t2.*t42.*t99;
t721 = t42.*t43.*t98;
t791 = t2.*t99.*t104;
t158 = dxi1.*t112;
t159 = dxi1.*t113;
t160 = dxi4.*t112;
t161 = dxi4.*t113;
t162 = dxi5.*t112;
t163 = dxi5.*t113;
t179 = dxi14.*t114;
t180 = dxi14.*t115;
t192 = mu3.*t107;
t193 = mu4.*t107;
t194 = mu3.*t108;
t195 = mu4.*t108;
t199 = mu4.*t111;
t216 = mu4.*t116;
t234 = dxi2.*t204;
t238 = dxi2.*t205;
t239 = dxi3.*t204;
t241 = dxi6.*t200;
t243 = dxi2.*t206;
t244 = dxi2.*t207;
t245 = dxi3.*t205;
t247 = dxi6.*t201;
t248 = dxi6.*t202;
t249 = dxi7.*t200;
t250 = dxi2.*t208;
t251 = dxi3.*t206;
t252 = dxi3.*t207;
t255 = dxi6.*t203;
t256 = dxi7.*t201;
t257 = dxi7.*t202;
t258 = dxi3.*t208;
t262 = dxi7.*t203;
t287 = dxi6.*t213;
t288 = dxi10.*t209;
t290 = dxi6.*t214;
t291 = dxi7.*t213;
t292 = dxi10.*t210;
t293 = dxi10.*t211;
t294 = dxi6.*t215;
t295 = dxi7.*t214;
t297 = dxi10.*t212;
t298 = dxi11.*t211;
t299 = dxi7.*t215;
t302 = dxi11.*t212;
t306 = mu4.*t167;
t307 = mu4.*t171;
t308 = mu4.*t173;
t309 = mu4.*t177;
t310 = mu4.*t178;
t311 = mu4.*t105.*xi14;
t313 = mu4.*t106.*xi14;
t315 = mu4.*t181;
t316 = mu4.*t182;
t317 = dxi6.*t118;
t318 = dxi7.*t118;
t319 = dxi6.*t156;
t320 = dxi7.*t156;
t321 = t217.*xi14;
t322 = t218.*xi14;
t323 = t219.*xi14;
t324 = t220.*xi14;
t325 = t221.*xi14;
t326 = t222.*xi14;
t327 = t223.*xi14;
t328 = t224.*xi14;
t329 = t225.*xi14;
t330 = t226.*xi14;
t335 = t118.*xi6;
t336 = t118.*xi7;
t337 = -t107;
t338 = -t108;
t339 = -t109;
t340 = -t110;
t341 = cos(t157);
t342 = sin(t157);
t349 = t116+xi10+xi11;
t386 = dxi1.*mu4.*t155;
t390 = dxi10.*mu3.*t155;
t393 = dxi10.*mu4.*t155;
t396 = dxi11.*mu4.*t155;
t426 = dxi14.*mu4.*t156;
t435 = t32.*t117;
t436 = t33.*t117;
t437 = t34.*t117;
t438 = mu4.*t118.*xi14;
t441 = dxi1.*t156.*xi6;
t442 = dxi1.*t156.*xi7;
t443 = dxi4.*t156.*xi6;
t444 = dxi4.*t156.*xi7;
t445 = dxi5.*t156.*xi6;
t446 = dxi5.*t156.*xi7;
t447 = -t187;
t448 = -t188;
t449 = -t189;
t450 = -t190;
t451 = -t191;
t452 = -t196;
t453 = -t197;
t454 = -t198;
t458 = t58.*t117;
t459 = t63.*t117;
t461 = dxi10.*t23.*t117;
t462 = t68.*t117;
t465 = dxi10.*t24.*t117;
t466 = dxi10.*t25.*t117;
t469 = t73.*t117;
t470 = dxi10.*t26.*t117;
t471 = dxi11.*t25.*t117;
t473 = t77.*t117;
t474 = t78.*t117;
t475 = dxi11.*t26.*t117;
t476 = t80.*t117;
t477 = t81.*t117;
t478 = t82.*t117;
t479 = t83.*t117;
t480 = t84.*t117;
t481 = t85.*t117;
t482 = t86.*t117;
t483 = t87.*t117;
t484 = t89.*t117;
t485 = mu4.*t6.*t118;
t486 = mu4.*t10.*t118;
t487 = mu4.*t12.*t118;
t489 = dxi14.*t30.*t118;
t491 = dxi14.*t31.*t118;
t492 = mu4.*t16.*t118;
t493 = mu4.*t17.*t118;
t494 = mu4.*t18.*t118;
t495 = mu4.*t19.*t118;
t498 = t95.*t121;
t499 = t95.*t122;
t500 = t95.*t123;
t501 = t95.*t124;
t502 = t95.*t125;
t503 = t95.*t126;
t504 = t95.*t127;
t505 = t95.*t128;
t506 = t95.*t129;
t507 = t95.*t130;
t508 = t95.*t131;
t509 = t95.*t132;
t510 = t95.*t134;
t511 = t95.*t135;
t512 = t95.*t136;
t513 = t95.*t138;
t514 = t95.*t139;
C_xi = ft_1({dxi1,dxi10,dxi11,dxi12,dxi13,dxi14,dxi2,dxi3,dxi4,dxi5,mu1,mu2,mu3,mu4,t10,t100,t101,t104,t105,t106,t109,t11,t110,t111,t112,t113,t114,t115,t117,t118,t119,t12,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t13,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t14,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t15,t150,t151,t152,t153,t154,t155,t156,t158,t159,t16,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t17,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t18,t180,t181,t182,t183,t184,t185,t186,t19,t192,t193,t194,t195,t199,t2,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t21,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t22,t220,t221,t222,t223,t224,t225,t226,t227,t23,t232,t233,t234,t235,t236,t237,t238,t239,t24,t240,t241,t242,t243,t244,t245,t247,t248,t249,t25,t250,t251,t252,t255,t256,t257,t258,t26,t262,t27,t28,t287,t288,t29,t290,t291,t292,t293,t294,t295,t297,t298,t299,t3,t30,t302,t306,t307,t308,t309,t31,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t32,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t33,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t34,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t38,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t39,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t4,t40,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t41,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t42,t426,t431,t432,t433,t434,t435,t436,t437,t438,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t456,t457,t458,t459,t461,t462,t465,t466,t470,t471,t475,t485,t486,t487,t489,t491,t492,t493,t494,t495,t498,t499,t5,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t6,t62,t67,t7,t70,t71,t716,t717,t718,t719,t72,t721,t74,t75,t76,t79,t791,t8,t88,t9,t90,t91,t92,t93,t94,t95,t96,t97,t98,xi14,xi2,xi3,xi6,xi7});
end
function C_xi = ft_1(ct)
dxi1 = ct{1};
dxi10 = ct{2};
dxi11 = ct{3};
dxi12 = ct{4};
dxi13 = ct{5};
dxi14 = ct{6};
dxi2 = ct{7};
dxi3 = ct{8};
dxi4 = ct{9};
dxi5 = ct{10};
mu1 = ct{11};
mu2 = ct{12};
mu3 = ct{13};
mu4 = ct{14};
t10 = ct{15};
t100 = ct{16};
t101 = ct{17};
t104 = ct{18};
t105 = ct{19};
t106 = ct{20};
t109 = ct{21};
t11 = ct{22};
t110 = ct{23};
t111 = ct{24};
t112 = ct{25};
t113 = ct{26};
t114 = ct{27};
t115 = ct{28};
t117 = ct{29};
t118 = ct{30};
t119 = ct{31};
t12 = ct{32};
t120 = ct{33};
t121 = ct{34};
t122 = ct{35};
t123 = ct{36};
t124 = ct{37};
t125 = ct{38};
t126 = ct{39};
t127 = ct{40};
t128 = ct{41};
t129 = ct{42};
t13 = ct{43};
t130 = ct{44};
t131 = ct{45};
t132 = ct{46};
t133 = ct{47};
t134 = ct{48};
t135 = ct{49};
t136 = ct{50};
t137 = ct{51};
t138 = ct{52};
t139 = ct{53};
t14 = ct{54};
t140 = ct{55};
t141 = ct{56};
t142 = ct{57};
t143 = ct{58};
t144 = ct{59};
t145 = ct{60};
t146 = ct{61};
t147 = ct{62};
t148 = ct{63};
t149 = ct{64};
t15 = ct{65};
t150 = ct{66};
t151 = ct{67};
t152 = ct{68};
t153 = ct{69};
t154 = ct{70};
t155 = ct{71};
t156 = ct{72};
t158 = ct{73};
t159 = ct{74};
t16 = ct{75};
t160 = ct{76};
t161 = ct{77};
t162 = ct{78};
t163 = ct{79};
t164 = ct{80};
t165 = ct{81};
t166 = ct{82};
t167 = ct{83};
t168 = ct{84};
t169 = ct{85};
t17 = ct{86};
t170 = ct{87};
t171 = ct{88};
t172 = ct{89};
t173 = ct{90};
t174 = ct{91};
t175 = ct{92};
t176 = ct{93};
t177 = ct{94};
t178 = ct{95};
t179 = ct{96};
t18 = ct{97};
t180 = ct{98};
t181 = ct{99};
t182 = ct{100};
t183 = ct{101};
t184 = ct{102};
t185 = ct{103};
t186 = ct{104};
t19 = ct{105};
t192 = ct{106};
t193 = ct{107};
t194 = ct{108};
t195 = ct{109};
t199 = ct{110};
t2 = ct{111};
t200 = ct{112};
t201 = ct{113};
t202 = ct{114};
t203 = ct{115};
t204 = ct{116};
t205 = ct{117};
t206 = ct{118};
t207 = ct{119};
t208 = ct{120};
t209 = ct{121};
t21 = ct{122};
t210 = ct{123};
t211 = ct{124};
t212 = ct{125};
t213 = ct{126};
t214 = ct{127};
t215 = ct{128};
t216 = ct{129};
t217 = ct{130};
t218 = ct{131};
t219 = ct{132};
t22 = ct{133};
t220 = ct{134};
t221 = ct{135};
t222 = ct{136};
t223 = ct{137};
t224 = ct{138};
t225 = ct{139};
t226 = ct{140};
t227 = ct{141};
t23 = ct{142};
t232 = ct{143};
t233 = ct{144};
t234 = ct{145};
t235 = ct{146};
t236 = ct{147};
t237 = ct{148};
t238 = ct{149};
t239 = ct{150};
t24 = ct{151};
t240 = ct{152};
t241 = ct{153};
t242 = ct{154};
t243 = ct{155};
t244 = ct{156};
t245 = ct{157};
t247 = ct{158};
t248 = ct{159};
t249 = ct{160};
t25 = ct{161};
t250 = ct{162};
t251 = ct{163};
t252 = ct{164};
t255 = ct{165};
t256 = ct{166};
t257 = ct{167};
t258 = ct{168};
t26 = ct{169};
t262 = ct{170};
t27 = ct{171};
t28 = ct{172};
t287 = ct{173};
t288 = ct{174};
t29 = ct{175};
t290 = ct{176};
t291 = ct{177};
t292 = ct{178};
t293 = ct{179};
t294 = ct{180};
t295 = ct{181};
t297 = ct{182};
t298 = ct{183};
t299 = ct{184};
t3 = ct{185};
t30 = ct{186};
t302 = ct{187};
t306 = ct{188};
t307 = ct{189};
t308 = ct{190};
t309 = ct{191};
t31 = ct{192};
t310 = ct{193};
t311 = ct{194};
t312 = ct{195};
t313 = ct{196};
t314 = ct{197};
t315 = ct{198};
t316 = ct{199};
t317 = ct{200};
t318 = ct{201};
t319 = ct{202};
t32 = ct{203};
t320 = ct{204};
t321 = ct{205};
t322 = ct{206};
t323 = ct{207};
t324 = ct{208};
t325 = ct{209};
t326 = ct{210};
t327 = ct{211};
t328 = ct{212};
t329 = ct{213};
t33 = ct{214};
t330 = ct{215};
t331 = ct{216};
t332 = ct{217};
t333 = ct{218};
t334 = ct{219};
t335 = ct{220};
t336 = ct{221};
t337 = ct{222};
t338 = ct{223};
t339 = ct{224};
t34 = ct{225};
t340 = ct{226};
t341 = ct{227};
t342 = ct{228};
t343 = ct{229};
t344 = ct{230};
t345 = ct{231};
t346 = ct{232};
t347 = ct{233};
t348 = ct{234};
t349 = ct{235};
t356 = ct{236};
t357 = ct{237};
t358 = ct{238};
t359 = ct{239};
t360 = ct{240};
t361 = ct{241};
t362 = ct{242};
t363 = ct{243};
t364 = ct{244};
t365 = ct{245};
t366 = ct{246};
t367 = ct{247};
t368 = ct{248};
t369 = ct{249};
t370 = ct{250};
t371 = ct{251};
t372 = ct{252};
t373 = ct{253};
t374 = ct{254};
t375 = ct{255};
t376 = ct{256};
t377 = ct{257};
t378 = ct{258};
t379 = ct{259};
t38 = ct{260};
t380 = ct{261};
t381 = ct{262};
t382 = ct{263};
t383 = ct{264};
t384 = ct{265};
t385 = ct{266};
t386 = ct{267};
t387 = ct{268};
t388 = ct{269};
t389 = ct{270};
t39 = ct{271};
t390 = ct{272};
t391 = ct{273};
t392 = ct{274};
t393 = ct{275};
t394 = ct{276};
t395 = ct{277};
t396 = ct{278};
t397 = ct{279};
t398 = ct{280};
t399 = ct{281};
t4 = ct{282};
t40 = ct{283};
t400 = ct{284};
t401 = ct{285};
t402 = ct{286};
t403 = ct{287};
t404 = ct{288};
t405 = ct{289};
t406 = ct{290};
t407 = ct{291};
t408 = ct{292};
t409 = ct{293};
t41 = ct{294};
t410 = ct{295};
t411 = ct{296};
t412 = ct{297};
t413 = ct{298};
t414 = ct{299};
t415 = ct{300};
t416 = ct{301};
t417 = ct{302};
t418 = ct{303};
t419 = ct{304};
t42 = ct{305};
t426 = ct{306};
t431 = ct{307};
t432 = ct{308};
t433 = ct{309};
t434 = ct{310};
t435 = ct{311};
t436 = ct{312};
t437 = ct{313};
t438 = ct{314};
t441 = ct{315};
t442 = ct{316};
t443 = ct{317};
t444 = ct{318};
t445 = ct{319};
t446 = ct{320};
t447 = ct{321};
t448 = ct{322};
t449 = ct{323};
t450 = ct{324};
t451 = ct{325};
t452 = ct{326};
t453 = ct{327};
t454 = ct{328};
t456 = ct{329};
t457 = ct{330};
t458 = ct{331};
t459 = ct{332};
t461 = ct{333};
t462 = ct{334};
t465 = ct{335};
t466 = ct{336};
t470 = ct{337};
t471 = ct{338};
t475 = ct{339};
t485 = ct{340};
t486 = ct{341};
t487 = ct{342};
t489 = ct{343};
t491 = ct{344};
t492 = ct{345};
t493 = ct{346};
t494 = ct{347};
t495 = ct{348};
t498 = ct{349};
t499 = ct{350};
t5 = ct{351};
t500 = ct{352};
t501 = ct{353};
t502 = ct{354};
t503 = ct{355};
t504 = ct{356};
t505 = ct{357};
t506 = ct{358};
t507 = ct{359};
t508 = ct{360};
t509 = ct{361};
t510 = ct{362};
t511 = ct{363};
t512 = ct{364};
t513 = ct{365};
t514 = ct{366};
t6 = ct{367};
t62 = ct{368};
t67 = ct{369};
t7 = ct{370};
t70 = ct{371};
t71 = ct{372};
t716 = ct{373};
t717 = ct{374};
t718 = ct{375};
t719 = ct{376};
t72 = ct{377};
t721 = ct{378};
t74 = ct{379};
t75 = ct{380};
t76 = ct{381};
t79 = ct{382};
t791 = ct{383};
t8 = ct{384};
t88 = ct{385};
t9 = ct{386};
t90 = ct{387};
t91 = ct{388};
t92 = ct{389};
t93 = ct{390};
t94 = ct{391};
t95 = ct{392};
t96 = ct{393};
t97 = ct{394};
t98 = ct{395};
xi14 = ct{396};
xi2 = ct{397};
xi3 = ct{398};
xi6 = ct{399};
xi7 = ct{400};
t515 = t95.*t140;
t516 = t95.*t142;
t517 = t96.*t126;
t518 = t96.*t127;
t519 = t96.*t128;
t520 = t96.*t129;
t521 = t96.*t131;
t522 = t96.*t133;
t523 = t96.*t134;
t524 = t96.*t135;
t525 = t96.*t136;
t526 = t96.*t137;
t527 = t96.*t138;
t528 = t96.*t139;
t529 = t96.*t140;
t530 = t96.*t141;
t531 = t96.*t142;
t532 = t96.*t143;
t533 = t96.*t144;
t534 = t96.*t145;
t535 = t96.*t146;
t536 = t96.*t147;
t537 = t96.*t148;
t538 = t96.*t149;
t539 = t96.*t150;
t540 = t96.*t151;
t541 = t96.*t152;
t542 = t96.*t153;
t543 = t96.*t154;
t553 = mu4.*t227.*xi14;
t554 = dxi1.*t2.*t155;
t555 = t3.*t23.*t155;
t556 = t3.*t24.*t155;
t557 = t3.*t25.*t155;
t558 = t4.*t25.*t155;
t559 = t3.*t26.*t155;
t560 = t4.*t26.*t155;
t561 = t5.*t23.*t155;
t562 = t5.*t24.*t155;
t563 = t5.*t25.*t155;
t564 = t8.*t23.*t155;
t565 = t7.*t25.*t155;
t566 = t5.*t26.*t155;
t567 = t8.*t24.*t155;
t568 = t8.*t25.*t155;
t569 = t7.*t26.*t155;
t570 = t9.*t25.*t155;
t571 = t8.*t26.*t155;
t572 = t9.*t26.*t155;
t573 = t11.*t23.*t155;
t574 = t11.*t24.*t155;
t575 = t11.*t25.*t155;
t576 = t14.*t23.*t155;
t577 = t13.*t25.*t155;
t578 = t11.*t26.*t155;
t579 = t14.*t24.*t155;
t580 = t14.*t25.*t155;
t581 = t13.*t26.*t155;
t582 = t15.*t25.*t155;
t583 = t14.*t26.*t155;
t584 = t15.*t26.*t155;
t586 = t6.*t30.*t156;
t587 = t6.*t31.*t156;
t589 = t10.*t30.*t156;
t590 = t10.*t31.*t156;
t591 = t12.*t30.*t156;
t592 = t12.*t31.*t156;
t593 = t16.*t30.*t156;
t594 = t16.*t31.*t156;
t595 = t17.*t30.*t156;
t596 = t17.*t31.*t156;
t597 = t18.*t30.*t156;
t598 = t18.*t31.*t156;
t599 = t19.*t30.*t156;
t600 = t19.*t31.*t156;
t604 = t98.*t119.*xi6;
t605 = t98.*t120.*xi6;
t606 = t98.*t121.*xi6;
t607 = t98.*t121.*xi7;
t608 = t98.*t122.*xi6;
t609 = t98.*t123.*xi6;
t610 = t98.*t122.*xi7;
t611 = t98.*t123.*xi7;
t612 = t98.*t124.*xi6;
t613 = -t365;
t614 = t98.*t124.*xi7;
t615 = -t367;
t616 = -t368;
t617 = -t369;
t618 = -t370;
t619 = -t371;
t620 = -t372;
t621 = -t373;
t622 = -t374;
t623 = -t375;
t624 = -t376;
t625 = -t377;
t626 = -t378;
t627 = -t379;
t628 = -t380;
t629 = -t381;
t630 = -t382;
t631 = -t383;
t632 = -t384;
t633 = -t385;
t634 = -t387;
t635 = -t388;
t636 = -t389;
t638 = -t391;
t639 = -t392;
t641 = -t394;
t642 = -t395;
t644 = -t397;
t645 = -t398;
t646 = -t399;
t647 = -t400;
t648 = -t401;
t649 = -t402;
t650 = -t403;
t651 = -t404;
t652 = -t405;
t653 = -t406;
t654 = -t407;
t655 = -t408;
t656 = -t409;
t657 = -t410;
t658 = -t411;
t659 = -t412;
t660 = -t413;
t661 = -t414;
t662 = -t415;
t663 = -t416;
t664 = -t417;
t665 = -t418;
t666 = -t419;
t667 = t100.*t137.*xi14;
t668 = t100.*t141.*xi14;
t669 = t100.*t144.*xi14;
t670 = t100.*t146.*xi14;
t671 = t100.*t147.*xi14;
t672 = t100.*t148.*xi14;
t674 = t100.*t150.*xi14;
t675 = t100.*t152.*xi14;
t676 = t100.*t153.*xi14;
t677 = t100.*t154.*xi14;
t678 = -t431;
t679 = -t432;
t680 = -t433;
t681 = -t434;
t682 = mu4.*t42.*t156;
t694 = t117.*t133;
t695 = t117.*t137;
t696 = t117.*t141;
t697 = t117.*t143;
t698 = t117.*t144;
t699 = t117.*t145;
t700 = t117.*t146;
t701 = t117.*t147;
t702 = t117.*t148;
t703 = t117.*t149;
t704 = t117.*t150;
t705 = t117.*t151;
t706 = t117.*t152;
t707 = t117.*t153;
t708 = t117.*t154;
t790 = mu4.*t100.*t227;
t792 = mu4.*t104.*t155;
t793 = -t716;
t794 = -t718;
t795 = t114+t115+xi14;
t796 = -t721;
t811 = mu4.*t156.*t227;
t812 = t2.*t104.*t155;
t814 = -t791;
t817 = t39+t40+t345;
t822 = t346+t347;
t825 = t41+t343+t344;
t848 = t111+t167+t171+t173+t177+t178+t181+t182;
t852 = dxi14+t105+t106+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226;
t350 = dxi1.*t335;
t351 = dxi1.*t336;
t352 = dxi4.*t335;
t353 = dxi4.*t336;
t354 = dxi5.*t335;
t355 = dxi5.*t336;
t439 = -t179;
t440 = -t180;
t455 = -t199;
t460 = dxi2.*t435;
t463 = dxi2.*t436;
t464 = dxi3.*t435;
t467 = dxi2.*t437;
t468 = dxi3.*t436;
t472 = dxi3.*t437;
t488 = mu4.*t317.*xi14;
t490 = mu4.*t318.*xi14;
t496 = dxi2.*t341;
t497 = dxi3.*t341;
t544 = -t306;
t545 = -t307;
t546 = -t308;
t547 = -t309;
t548 = -t310;
t549 = -t315;
t550 = -t316;
t551 = -t319;
t552 = -t320;
t585 = dxi1.*mu4.*t342;
t588 = dxi14.*mu4.*t342;
t601 = mu4.*t341.*xi14;
t602 = dxi1.*t342.*xi2;
t603 = dxi1.*t342.*xi3;
t637 = -t390;
t640 = -t393;
t643 = -t396;
t673 = -t426;
t683 = mu4.*t6.*t341;
t685 = dxi14.*t25.*t341;
t687 = dxi14.*t26.*t341;
t688 = mu4.*t10.*t341;
t689 = mu4.*t12.*t341;
t690 = mu4.*t16.*t341;
t691 = mu4.*t17.*t341;
t692 = mu4.*t18.*t341;
t693 = mu4.*t19.*t341;
t709 = -t485;
t710 = -t486;
t711 = -t487;
t712 = -t492;
t713 = -t493;
t714 = -t494;
t715 = -t495;
t720 = t104.*t216;
t722 = -t553;
t723 = t42.*t438;
t724 = t6.*t25.*t342;
t725 = t6.*t26.*t342;
t726 = t10.*t25.*t342;
t727 = t10.*t26.*t342;
t728 = t12.*t25.*t342;
t729 = t12.*t26.*t342;
t730 = t16.*t25.*t342;
t731 = t16.*t26.*t342;
t732 = t17.*t25.*t342;
t733 = t17.*t26.*t342;
t734 = t18.*t25.*t342;
t735 = t18.*t26.*t342;
t736 = t19.*t25.*t342;
t737 = t19.*t26.*t342;
t738 = -t555;
t739 = -t556;
t740 = -t557;
t741 = -t558;
t742 = -t559;
t743 = -t560;
t744 = -t561;
t745 = -t562;
t746 = -t563;
t747 = -t564;
t748 = -t565;
t749 = -t566;
t750 = -t567;
t751 = -t568;
t752 = -t569;
t753 = -t570;
t754 = -t571;
t755 = -t572;
t756 = -t573;
t757 = -t574;
t758 = -t575;
t759 = -t576;
t760 = -t577;
t761 = -t578;
t762 = -t579;
t763 = -t580;
t764 = -t581;
t765 = -t582;
t766 = -t583;
t767 = -t584;
t768 = -t586;
t769 = -t587;
t771 = -t589;
t772 = -t590;
t773 = -t591;
t774 = -t592;
t775 = -t593;
t776 = -t594;
t777 = -t595;
t778 = -t596;
t779 = -t597;
t780 = -t598;
t781 = -t599;
t782 = -t600;
t813 = -t790;
t815 = -t792;
t816 = t32+t33+t34+t216;
t818 = mu4.*t104.*t349;
t819 = mu4.*t227.*t342;
t820 = -t811;
t821 = -t812;
t824 = t38+t817;
t826 = t112+t113+t349;
t827 = (dxi11.*t817)./2.0;
t828 = (dxi12.*t822)./2.0;
t829 = (dxi13.*t822)./2.0;
t830 = mu4.*t227.*t795;
t834 = (dxi14.*t825)./2.0;
t837 = t335+t336+t795;
t838 = t213+t214+t215+t438;
t849 = mu4.*t848;
t853 = t3+t4+t5+t7+t8+t9+t11+t13+t14+t15+t848;
t854 = mu4.*t852.*xi14;
t858 = t6+t10+t12+t16+t17+t18+t19+t164+t165+t166+t168+t169+t170+t172+t174+t175+t176+t339+t340;
t869 = t317+t318+t441+t442+t443+t444+t445+t446+t852;
t684 = mu4.*t496.*xi14;
t686 = mu4.*t497.*xi14;
t770 = -t588;
t783 = -t683;
t784 = -t688;
t785 = -t689;
t786 = -t690;
t787 = -t691;
t788 = -t692;
t789 = -t693;
t797 = -t724;
t798 = -t725;
t799 = -t726;
t800 = -t727;
t801 = -t728;
t802 = -t729;
t803 = -t730;
t804 = -t731;
t805 = -t732;
t806 = -t733;
t807 = -t734;
t808 = -t735;
t809 = -t736;
t810 = -t737;
t823 = -t819;
t831 = t104.*t816;
t832 = -t828;
t833 = -t829;
t835 = (dxi10.*t824)./2.0;
t836 = mu4.*t104.*t826;
t839 = t42.*t838;
t840 = mu4.*t227.*t837;
t841 = t28+t29+t30+t31+t838;
t842 = t458+t459+t462+t683;
t843 = t209+t210+t211+t212+t816;
t850 = -t849;
t855 = mu4.*t853;
t857 = t673+t709+t710+t711+t712+t713+t714+t715;
t859 = mu4.*t858;
t863 = t101+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t439+t440;
t866 = t158+t159+t160+t161+t162+t163+t337+t338+t853;
t870 = mu4.*t869.*xi14;
t871 = t133+t137+t141+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t455+t544+t545+t546+t547+t548+t549+t550;
t872 = t350+t351+t352+t353+t354+t355+t551+t552+t858;
t875 = t496+t497+t602+t603+t869;
t882 = t94+t312+t314+t489+t491+t667+t668+t669+t670+t671+t672+t674+t675+t676+t677+t678+t679+t680+t681+t768+t769+t771+t772+t773+t774+t775+t776+t777+t778+t779+t780+t781+t782;
t883 = t88+t90+t91+t92+t93+t94+t287+t290+t291+t294+t295+t299+t311+t312+t313+t314+t387+t388+t389+t391+t392+t394+t395+t397+t398+t399+t400+t401+t402+t403+t404+t405+t406+t407+t488+t490+t586+t587+t589+t590+t591+t592+t678+t679+t680+t681;
t884 = t62+t67+t70+t71+t72+t74+t75+t76+t79+t88+t90+t91+t92+t93+t94+t287+t288+t290+t291+t292+t293+t294+t295+t297+t298+t299+t302+t311+t312+t313+t314+t488+t489+t490+t491+t655+t656+t657+t658+t659+t660+t661+t662+t663+t664+t665+t666+t678+t679+t680+t681+t775+t776+t777+t778+t779+t780+t781+t782;
t888 = t88+t90+t91+t92+t93+t94+t288+t292+t293+t297+t298+t302+t311+t312+t313+t314+t489+t491+t634+t635+t636+t638+t639+t641+t642+t644+t645+t646+t647+t648+t649+t650+t651+t652+t653+t654+t655+t656+t657+t658+t659+t660+t661+t662+t663+t664+t665+t666+t678+t679+t680+t681+t768+t769+t771+t772+t773+t774+t775+t776+t777+t778+t779+t780+t781+t782;
t844 = t27+t841;
t845 = t42.*t841;
t847 = t104.*t843;
t851 = t232+t233+t236+t237+t242+t842;
t856 = -t855;
t860 = -t859;
t861 = t770+t783+t784+t785+t786+t787+t788+t789;
t862 = t827+t832+t833+t834+t835;
t864 = mu4.*t863;
t867 = mu4.*t866;
t873 = mu4.*t872;
t876 = mu4.*t875.*xi14;
t877 = t452+t453+t454+t522+t526+t530+t532+t533+t534+t535+t536+t537+t538+t539+t540+t541+t542+t543+t857;
t878 = t192+t193+t194+t195+t517+t518+t519+t520+t521+t523+t524+t525+t527+t528+t529+t531+t871;
t886 = t460+t463+t464+t467+t468+t472+t555+t556+t557+t558+t559+t560+t684+t686+t724+t725+t883;
t887 = t685+t687+t797+t798+t799+t800+t801+t802+t803+t804+t805+t806+t807+t808+t809+t810+t882;
t889 = t234+t238+t239+t243+t244+t245+t250+t251+t252+t258+t356+t357+t358+t359+t360+t361+t362+t363+t364+t366+t460+t463+t464+t467+t468+t472+t555+t556+t557+t558+t559+t560+t684+t686+t724+t725+t884;
t890 = t461+t465+t466+t470+t471+t475+t685+t687+t738+t739+t740+t741+t742+t743+t744+t745+t746+t747+t748+t749+t750+t751+t752+t753+t754+t755+t756+t757+t758+t759+t760+t761+t762+t763+t764+t765+t766+t767+t797+t798+t799+t800+t801+t802+t803+t804+t805+t806+t807+t808+t809+t810+t888;
t891 = t235+t240+t241+t247+t248+t249+t255+t256+t257+t262+t461+t465+t466+t470+t471+t475+t604+t605+t606+t607+t608+t609+t610+t611+t612+t613+t614+t615+t616+t617+t618+t619+t620+t621+t622+t623+t624+t625+t626+t627+t628+t629+t630+t631+t632+t633+t685+t687+t738+t739+t740+t741+t742+t743+t744+t745+t746+t747+t748+t749+t750+t751+t752+t753+t754+t755+t756+t757+t758+t759+t760+t761+t762+t763+t764+t765+t766+t767+t797+t798+t799+t800+t801+t802+t803+t804+t805+t806+t807+t808+t809+t810+t884;
t846 = t42.*t844;
t865 = -t864;
t868 = -t867;
t874 = -t873;
t879 = t637+t640+t643+t694+t695+t696+t697+t698+t699+t700+t701+t702+t703+t704+t705+t706+t707+t708+t861;
t880 = t126+t127+t128+t129+t131+t134+t135+t136+t138+t139+t140+t142+t877;
t881 = t125+t130+t132+t880;
t885 = t447+t448+t449+t450+t451+t502+t503+t504+t505+t506+t507+t508+t509+t510+t511+t512+t513+t514+t515+t516+t879;
mt1 = [t234+t235+t238+t239+t240+t241+t243+t244+t245+t247+t248+t249+t250+t251+t252+t255+t256+t257+t258+t262+t460+t461+t463+t464+t465+t466+t467+t468+t470+t471+t472+t475+t613+t615+t616+t617+t618+t619+t620+t621+t622+t623+t624+t625+t626+t627+t628+t629+t630+t631+t632+t633+t684+t685+t686+t687+t744+t745+t746+t747+t748+t749+t750+t751+t752+t753+t754+t755+t756+t757+t758+t759+t760+t761+t762+t763+t764+t765+t766+t767+t799+t800+t801+t802+t803+t804+t805+t806+t807+t808+t809+t810+t884+dxi2.*t21+dxi2.*t22+dxi3.*t21+dxi2.*t23+dxi3.*t22+dxi2.*t24+dxi3.*t23+dxi2.*t25+dxi3.*t24+dxi2.*t26+dxi3.*t25+dxi3.*t26+dxi2.*mu1.*xi2,t119+t120+t121+t122+t123+t124+t885-dxi1.*mu1.*xi2,t119+t120+t121+t122+t123+t124+t885,t889,t889,t183+t184+t185+t186+t498+t499+t500+t501+t881+t95.*t119+t95.*t120+dxi2.*mu2.*t98+dxi3.*mu2.*t98,t183+t184+t185+t186+t498+t499+t500+t501+t880,t886,t886];
mt2 = [t878+t117.*t121+t117.*t122+t117.*t123+t117.*t124+dxi2.*mu3.*t155+dxi2.*mu4.*t155+dxi3.*mu3.*t155+dxi3.*mu4.*t155,-mu4.*(t866-dxi2.*t155-dxi3.*t155+dxi1.*t117.*xi2+dxi1.*t117.*xi3),t876,t876,-mu4.*t6-mu4.*t10-mu4.*t12-mu4.*t16-mu4.*t17-mu4.*t18-mu4.*t19+mu4.*t109+mu4.*t110+mu4.*t319+mu4.*t320+t97.*t137+t97.*t141+t97.*t144+t97.*t146+t97.*t147+t97.*t148+t118.*t128+t97.*t150+t118.*t129+t97.*t152+t97.*t153+t97.*t154+t118.*t135+t118.*t138+t118.*t140+t118.*t142+t123.*t341+t124.*t341+dxi2.*mu4.*t342+dxi3.*mu4.*t342,0.0,0.0,dxi1.*(t21+t22+t23+t24+t25+t26+t204+t205+t206+t207+t208+t435+t436+t437+t601+mu1.*xi2),0.0,0.0,t851,t851,t457,t348,t842,t842,t554,t386,t683];
mt3 = [t683,t585,0.0,0.0,dxi1.*(t21+t22+t23+t24+t25+t26+t204+t205+t206+t207+t208+t435+t436+t437+t601),0.0,0.0,t851,t851,t457,t348,t842,t842,t554,t386,t683,t683,t585,0.0,0.0,t891,t885,t885,t884,t884,t881,t880,t883,t883,t878,t868,t870,t870,t874,0.0,0.0,t891,t885,t885,t884,t884,t881,t880,t883,t883,t878,t868,t870,t870,t874,0.0,0.0,t42.*(t200+t201+t202+t203+t844+t21.*t95+t22.*t95),t796,t796,t846,t846,0.0,0.0,t839,t839,t719,t456,t723,t723,t682,0.0,0.0,t42.*(t200+t201+t202+t203+t841),t794,t794,t845,t845,0.0,0.0,t839,t839,t719,t456,t723,t723,t682,0.0,0.0,t890,t879,t879,t888,t888,t877,t877,t862,t862,t871,t856,t854,t854,t860,0.0,0.0,t890,t879,t879,t888,t888,t877,t877,t862,t862,t871,t856,t854,t854,t860,0.0,0.0,t104.*(t843+t23.*t117+t24.*t117+t25.*t117+t26.*t117),t821,t821,t847,t847,t814,t814,t831,t831,0.0,0.0,t720,t720,t717,0.0,0.0,mu4.*t104.*(t826+t117.*xi2+t117.*xi3),t815,t815,t836,t836,t793];
mt4 = [t793,t818,t818,0.0,0.0,t720,t720,t717,0.0,0.0,t887,t861,t861,t882,t882,t857,t857,t865,t865,t850,t850,t94,t94,t722,0.0,0.0,t887,t861,t861,t882,t882,t857,t857,t865,t865,t850,t850,t94,t94,t722,0.0,0.0,mu4.*t227.*(t837+t341.*xi2+t341.*xi3),t823,t823,t840,t840,t820,t820,t830,t830,t813,t813,t553,t553,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
C_xi = reshape([mt1,mt2,mt3,mt4],16,16);
end
