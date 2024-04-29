function tau = get_tau_F(in1,in2)
%GET_TAU_F
%    TAU = GET_TAU_F(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    10-Apr-2024 16:26:21

F1 = in2(1,:);
F2 = in2(2,:);
F3 = in2(3,:);
F4 = in2(4,:);
F5 = in2(5,:);
F6 = in2(6,:);
q1 = in1(1,:);
q2 = in1(2,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q1);
t5 = sin(q2);
t6 = q1+q2;
t17 = sqrt(2.0);
t18 = atan(3.0./8.0);
t19 = atan(8.0./3.0);
t20 = pi./2.0;
t34 = sqrt(7.3e+1);
t7 = cos(t6);
t8 = sin(t6);
t9 = t2.*9.0;
t10 = t3.*9.0;
t11 = t3.*2.4e+1;
t12 = t3.*6.0e+1;
t13 = t5.*9.0;
t14 = t4.*2.4e+1;
t15 = t5.*2.4e+1;
t16 = t4.*6.0e+1;
t21 = t2./2.0;
t22 = t2./5.0;
t25 = t4./2.0;
t26 = t4./5.0;
t28 = t5.*1.6e+2;
t29 = -t20;
t32 = q2+t18;
t33 = t4.*(3.0./4.0e+1);
t35 = q1+t20;
t41 = t6+t19;
t42 = -t18;
t43 = -t19;
t45 = t6+t20;
t23 = -t9;
t24 = -t10;
t27 = -t14;
t30 = t7.*2.4e+1;
t31 = t8.*9.0;
t36 = t7./5.0;
t37 = t8./5.0;
t38 = -t33;
t39 = q1+t29;
t40 = cos(t32);
t44 = cos(t35);
t46 = sin(t35);
t47 = cos(t41);
t48 = cos(t45);
t49 = sin(t45);
t51 = t6+t29;
t55 = q2+t42;
t56 = t6+t43;
t50 = cos(t39);
t52 = sin(t39);
t53 = cos(t51);
t54 = sin(t51);
t57 = cos(t55);
t58 = cos(t56);
t59 = t44.*(3.0./4.0e+1);
t60 = t46.*(3.0./4.0e+1);
t63 = t48.*(3.0./4.0e+1);
t64 = t49.*(3.0./4.0e+1);
t65 = t14+t23+4.1e+1;
t68 = t23+t27+4.1e+1;
t71 = (t34.*t40)./4.0e+1;
t75 = t34.*t47.*(3.0./8.0e+2);
t61 = t50.*(3.0./4.0e+1);
t62 = t52.*(3.0./4.0e+1);
t66 = t53.*(3.0./4.0e+1);
t67 = t54.*(3.0./4.0e+1);
t69 = 1.0./sqrt(t65);
t70 = 1.0./sqrt(t68);
t72 = t22+t59;
t74 = (t34.*t57)./4.0e+1;
t76 = -t75;
t77 = t34.*t58.*(3.0./8.0e+2);
t79 = t21+t36+t63;
t73 = t22+t61;
t78 = -t77;
t80 = t21+t36+t66;
t81 = t33+t74+t76+2.41e+2./8.0e+2;
t82 = t38+t71+t78+2.41e+2./8.0e+2;
t83 = 1.0./sqrt(t81);
t84 = 1.0./sqrt(t82);
tau = [F1.*(t17.*t70.*t72.*(t26+t60).*2.0e+1-t17.*t70.*t72.*(t26+t60-3.0./4.0e+1).*2.0e+1)+F2.*(t17.*t69.*t73.*(t26+t62).*2.0e+1-t17.*t69.*t73.*(t26+t62+3.0./4.0e+1).*2.0e+1)+F5.*(t79.*t84.*(t25+t37+t64)-t79.*t84.*(t25+t37+t64-3.0./4.0e+1))+F6.*(t80.*t83.*(t25+t37+t67)-t80.*t83.*(t25+t37+t67+3.0./4.0e+1));F4.*t17.*(t11+t13).*1.0./sqrt(t15+t24+4.1e+1).*(-1.0./8.0e+1)-(F6.*t17.*1.0./sqrt(t16-t34.*t47.*3.0+t34.*t57.*2.0e+1+2.41e+2).*(t12-t28+t30+t31))./8.0e+1+(F3.*t17.*(t11-t13).*1.0./sqrt(-t15+t24+4.1e+1))./8.0e+1+(F5.*t17.*1.0./sqrt(-t16+t34.*t40.*2.0e+1-t34.*t58.*3.0+2.41e+2).*(t12+t28+t30-t31))./8.0e+1];