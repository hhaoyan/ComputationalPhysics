%%%%%%%%%%%%%%%%%%%%
% call me!
%%%%%%%%%%%%%%%%%%%%
function clock3
  t = [0:1/3600*2:12];
  y = cost_function_t(t);

  plot(t,y);
  xlabel('time/hour');
  ylabel('cost function');

  mt = sortrows([t(:) y(:)],2);

                                % search minimum in +/-1 minute
  deltaT = 1/3600*60;
  final_result = [];
  for j=1:22
    x = fminbnd(@cost_function_t, mt(j,1)-deltaT, mt(j,1)+deltaT);
    final_result = [final_result;x cost_function_t(x)];
  end
  final_result = sortrows(final_result, 2);
  x1 = final_result(1,1); x2=final_result(2,1);
  y1 = final_result(1,2); y2=final_result(2,2);

  display('using fminbnd@MATLAB to find the actual solution as:');
  display(sprintf('%.10f %.10f', x1, x2));
  display('cost function at solution is:');
  display(sprintf('%.10f %.10f', y1, y2));

  function cost=cost_function_t(t)
    cost = cost_function(t/12*2*pi, t*2*pi, t*60*2*pi);
  end

  function cost = cost_function(h,m,s)
    h = mod(h,2*pi);
    m = mod(m,2*pi);
    s = mod(s,2*pi);

    sorted = sort([h(:) m(:) s(:)], 2);
    t1 = sorted(:,1); t2 = sorted(:,2); t3 = sorted(:,3);

    cost = (t2-t1-2*pi/3).^2 + (t3-t2-2*pi/3).^2 + (t1-t3+4*pi/3).^2;
  end
end
