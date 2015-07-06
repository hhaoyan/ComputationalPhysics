function compute()
  close all;
  
  data = dlmread('pion-correlation-function.dat', '', 1);

  symmetric_data = cell(250,1);
  for i=1:250
    symmetric_data{i} = symmetric_piece(data((i-1)*64+1:i*64,:));
  end

%%%%%%%%%%%%%%%%%%%%
% problem a)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem a...');
  
  data_real = zeros(250, 33);
  data_image = zeros(250, 33);
  for i=1:250
    data_real(i,:) = symmetric_data{i}(:,1);
    data_image(i,:) = symmetric_data{i}(:,2);
  end
  avg_real = mean(data_real); avg_image = mean(data_image);
  delta_real = sqrt(var(data_real)/250); delta_image = sqrt(var(data_image)/250);

  relative_delta = delta_real ./ avg_real;

  coefficients = polyfit([0:32], relative_delta, 1);
  a = coefficients(1);
  b = coefficients(2);

  figure;
  hold on;
  plot([0:32],relative_delta,'.');
  plot([0:32], b+a*[0:32]);
  xlabel('t');
  ylabel('\Delta C(t) / avg(C(t))');
  title('Problem a: relative error-t');
  text(0.5, 0.3, sprintf('relative error = %f + %f t', b, a), 'Units','normalized');

%%%%%%%%%%%%%%%%%%%%
% problem b) jackknife
%%%%%%%%%%%%%%%%%%%%
  display('entering problem b...');
  
  resampled_avgs = jackknife(@mean, data_real, 1);
  resampled_meff_t = log(resampled_avgs(:,1:32) ./ resampled_avgs(:,2:33));
  resampled_length = size(resampled_meff_t, 1);
  meff_t = mean(resampled_meff_t);
  delta_meff_t = sqrt((resampled_length-1)^2/resampled_length*var(resampled_meff_t));
  figure;
  plot([0:31], delta_meff_t, '.');
  xlabel('t');
  ylabel('\Delta m_{eff}(t)');
  title('Problem b: errors of effetive mass-t');
  figure;
  errorbar([0:31], meff_t, delta_meff_t, '.k');
  xlabel('t');
  ylabel('m_{eff}(t)');
  title('Problem b: effetive mass-t and their errors');

%%%%%%%%%%%%%%%%%%%%
% problem c)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem c...');
  
  [meff, dmeff, chi2, pvalue, tstart, tend] = chi2_best_fit(meff_t, delta_meff_t, 4);
  display(sprintf('fitting result: meff: %f, dmeff: %f, chi2: %f, pvalue: %f, tstart: %d, tend: %d', meff, dmeff, chi2, pvalue, tstart-1, tend-1));
  hold on;
  line([5 30], [meff meff]);
  a = axis; miny = a(3); maxy = a(4);
  line([tstart tstart], [miny+0.15*(maxy-miny) maxy-0.15*(maxy-miny)]);
  line([tend tend], [miny+0.15*(maxy-miny) maxy-0.15*(maxy-miny)]);

%%%%%%%%%%%%%%%%%%%%
% problem d)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem d...');
  
  resampled_avgs_2 = jackknife(@mean, data_real, 1);
  resampled_meff_t_2 = acosh((resampled_avgs_2(:,1:31) + (resampled_avgs_2(:,3:33))) ./ 2 ./ resampled_avgs_2(:,2:32));
  resampled_length_2 = size(resampled_meff_t_2, 1);
  meff_t_2 = mean(resampled_meff_t_2);
  delta_meff_t_2 = sqrt((resampled_length_2-1)^2/resampled_length_2*var(resampled_meff_t_2));
  figure;
  errorbar([1:31], meff_t_2, delta_meff_t_2, '.k');
  xlabel('t');
  ylabel('m_{eff}(t)');
  title('Problem d: effetive mass-t and their errors');

  [meff_2, dmeff_2, chi2_2, pvalue_2, tstart_2, tend_2] = chi2_best_fit(meff_t_2, delta_meff_t_2, 4);
  display(sprintf('fitting result: meff: %f, dmeff: %f, chi2: %f, pvalue: %f, tstart: %d, tend: %d', meff_2, dmeff_2, chi2_2, pvalue_2, tstart_2, tend_2));
  hold on;
  line([5 30], [meff meff]);
  a = axis; miny = a(3); maxy = a(4);
  line([tstart_2 tstart_2], [miny+0.05*(maxy-miny) miny+0.3*(maxy-miny)]);
  line([tend_2 tend_2], [miny+0.05*(maxy-miny) miny+0.3*(maxy-miny)]);

%%%%%%%%%%%%%%%%%%%%
% problem e)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem e...');
  
  [bootstat, bootsam] = bootstrp(1000, @cov, data_real);
  cov_mean = mean(bootstat);
  cov_low_bound = zeros(1, size(cov_mean,2));
  cov_top_bound = zeros(1, size(cov_mean,2));
  for j=1:size(cov_mean,2)
    sorted_cov = sortrows(bootstat(:,j),1);
    cov_low_bound(j) = sorted_cov(160);
    cov_top_bound(j) = sorted_cov(840);
  end
  N = size(data_real,2);
  cov_mean = reshape(cov_mean, [N,N]);
  cov_low_bound = reshape(cov_low_bound, [N N]);
  cov_top_bound = reshape(cov_top_bound, [N N]);
  cov_delta = (cov_top_bound - cov_low_bound)/2;
% uncomment this to draw cov matrix.
  figure;
  errorbar([1:N*N], reshape(cov_mean, [N*N 1]), reshape(cov_low_bound, [N*N 1]), reshape(cov_top_bound, [N*N 1]));
  title('problem e: scatter plot of cov matrix');
  figure;
  imagesc(cov_mean);
  title('problem e: image for cov matrix');

  [bootstat, bootsam] = bootstrp(1000, @corrcoef, data_real);
  corrcoef_mean = mean(bootstat);
  corrcoef_low_bound = zeros(1, size(corrcoef_mean,2));
  corrcoef_top_bound = zeros(1, size(corrcoef_mean,2));
  for j=1:size(corrcoef_mean,2)
    sorted_corrcoef = sortrows(bootstat(:,j),1);
    corrcoef_low_bound(j) = sorted_corrcoef(160);
    corrcoef_top_bound(j) = sorted_corrcoef(840);
  end
  N = size(data_real,2);
  corrcoef_mean = reshape(corrcoef_mean, [N,N]);
  corrcoef_low_bound = reshape(corrcoef_low_bound, [N N]);
  corrcoef_top_bound = reshape(corrcoef_top_bound, [N N]);
  corrcoef_delta = (corrcoef_top_bound - corrcoef_low_bound)/2;
                                % rho 3,4
  display(sprintf('rho34: %f, delta_rho34: %f\nrho35: %f, delta_rho35: %f', corrcoef_mean(4,5), corrcoef_delta(4,5), corrcoef_mean(4,6), corrcoef_delta(4,6)));

%%%%%%%%%%%%%%%%%%%%
% problem f)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem f...');
  
  cov_mean_inv = myinv(cov_mean);

%%%%%%%%%%%%%%%%%%%%
% problem g)
%%%%%%%%%%%%%%%%%%%%
  display('entering problem g...');
  
  best_fit = [1000 0.14 1e99 -1 -1];
  for i=1:33
    for j=3:33
      if i+j > 33
        break;
      end

      tstart = i;
      tend = i+j;
      foo = @(x) non_linear_fit(x(1), x(2), avg_real, cov_mean_inv, tstart, tend);
      x = fminsearch(foo, [1000 0.1410], optimset('TolX', 1e-8));
      chi2 = foo(x);
      if chi2 < best_fit(3)
        best_fit = [x(1) x(2) chi2 tstart tend];
      end
    end
  end
  display(sprintf('best fit for the try: tstart: %d, tend: %d, a0: %.6f, mpi: %.6f, chi2: %f', best_fit(4)-1, best_fit(5)-1, best_fit(1), best_fit(2), best_fit(3)));
  

  function [a0 mpi] = local_nonlinearfit(values)
    model_func = @(b,x) b(1) .* (exp(-b(2)*x)+exp(-b(2)*(64-x)));
    x = 14:17;
    beta0 = [731.15835 0.142031];
    [beta, R, J, CovB, MSE, ErrorModelInfo] = nlinfit(x,values(15:18),model_func,beta0);
    a0 = beta(1);mpi=beta(2);
  end

  fitvalues = [];
  for i=1:250
    [a0 mpi] = local_nonlinearfit(data_real(i,:));
    fitvalues = [fitvalues; a0 mpi];
  end
  fitvalue_mean = mean(fitvalues);
  fitvalue_delta = sqrt(var(fitvalues));
  fitvalue_corrcoef = corrcoef(fitvalues);
  display('multiple fits:');
  display(sprintf('mean: (%f,%f), delta: (%f,%f), corrcoef:', fitvalue_mean(1), fitvalue_mean(2), fitvalue_delta(1), fitvalue_delta(2)));
  fitvalue_corrcoef

end

function chi2 = non_linear_fit(a0, mpi, avg, cov_inv, tstart, tend)
  chi2 = 0;
  for i=tstart:tend
    for j=tstart:tend
      c1 = a0 * (exp(-mpi*i)+exp(-mpi*(64-i)));
      c2 = a0 * (exp(-mpi*j)+exp(-mpi*(64-j)));
      chi2 = chi2 + (c1-avg(i))*cov_inv(i,j)*(c2-avg(j));
    end
  end
end

function inversed = myinv(m)
  N = size(m,1);
  [q, r] = qr(m);
                   % use backward substitution to find the inverse of r
  invr = zeros(N,N);
  for i=1:N
    b = zeros(N,1);
    b(i) = 1;
    invr(:,i) = back_substitution(r, b, N);
  end

  inversed = invr*transpose(q);

  function col = back_substitution(A, b, N)
    col = zeros(N, 1);
    col(N) = b(N) / A(N,N);
    for j=N-1:-1:1
      col(j) = (b(j)-A(j,j+1:N)*col(j+1:N))/A(j,j);
    end
  end
end

function [result, dresult, chi_squared, pvalue, tstart, tend] = chi2_best_fit(input_seq, dinput_seq, pt_lbound)
  input_seq = input_seq(:);
  dinput_seq = dinput_seq(:);
  time_start = 1;
  time_end = size(input_seq,1);

  best_fit = [0 9999999 9999 -1 -1];
  for a=time_start:time_end
    for b=pt_lbound-1:50
      tstart = a;
      tend = tstart+b;
      if tend > time_end
        break
      end
      [mpi, dmpi, chi_2] = chi2_fit(input_seq(tstart:tend), dinput_seq(tstart:tend));
      if (dmpi<best_fit(2)&&chi_2<best_fit(3)) || chi_2 < best_fit(3)
        best_fit = [mpi dmpi chi_2 tstart tend];
      end
    end
  end
  result = best_fit(1);
  dresult = best_fit(2);
  chi_squared = best_fit(3);
  tstart = best_fit(4);
  tend = best_fit(5);
  pvalue = chi2cdf(chi_squared, tend-tstart);
end

function [x,dx,chi_2] = chi2_fit(input, dinput)
  input = input(:);
  dinput = dinput(:);
  N = size(input,1);

  A = 1./dinput;
  y = input ./ dinput;
  x = transpose(A)*y/(transpose(A)*A);
                                %chi_2 = sum((y-A.*x).^2);
  chi_2 = sum(((input-x)./dinput).^2);
  dx = sqrt(chi2inv(chi2cdf(9,1), N)/(transpose(A)*A));
end

function result = symmetric_piece(piece)
  result = zeros(33, 2);

  for i=0:32
    if i==0
      result(i+1,:) = piece(1,2:3);
    else
      result(i+1,:) = 0.5*(piece(i+1,2:3) + piece(64-i+1,2:3));
    end
  end
end
