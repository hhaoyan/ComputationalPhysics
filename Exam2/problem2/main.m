%%%%%%%%%%%%%%%%%%%%
% call me!
%%%%%%%%%%%%%%%%%%%%
function main()
  q2 = 0.90;

  [x, fval, exitflag] = fsolve(@tosolve, q2);
  display(sprintf('problem solved, flag=%d, x=%.6f, fval=%.12f', exitflag, x, fval));
end

function total = tosolve(q2)
  total = zeta00(q2)/(pi^(3/2)) - 1 - 0.25*q2;
end

function total = zeta00(q2)
  total = zeta00_part_1(q2) + zeta00_part_2(q2) + zeta00_part_3(q2) + zeta00_part_4(q2);
end

function total = zeta00_part_4(q2)
  function s = innersum(x)
    s = zeros(size(x,1), size(x,2));
    tdints = threed_integers();
    for j=2:size(tdints,1)
      n_squared = sum((tdints(j,1:3)).^2);
      s = s + exp(x*q2).*exp(-1./x*n_squared*(pi^2)) * tdints(j,4);
    end
    s = x.^(-1.5).*s;
  end
  fun = @(x) innersum(x);
  total = integral(fun, 0, 1)*sqrt(pi/4);
end

function total = zeta00_part_3(q2)
  fun = @(x) x.^(-1.5).*(exp(x.*q2)-1);
  total = integral(fun, 0, 1) * pi / 2;
end

function total = zeta00_part_2(q2)
  total = -pi;
end

function total = zeta00_part_1(q2)
  total = 0;

  tdints = threed_integers();
  for j=1:size(tdints, 1)
    n_squared = sum((tdints(j,1:3)).^2);
    total = total + exp(q2-n_squared) / (n_squared-q2) * tdints(j,4);
  end

  total = total /sqrt(4*pi);
end

function result = threed_integers()
  tbl = [0 0 0 1;...
         0 0 1 6;...
         0 1 1 12;...
         1 1 1 8;...
         0 0 2 6;...
         0 1 2 24;...
         1 1 2 24;...
         0 2 2 12;...
         0 0 3 6;...
         1 2 2 24;...
         0 1 3 24;...
         1 1 3 24;...
         2 2 2 8;...
         0 2 3 24;...
         1 2 3 48;];

  result = tbl;
end
