function B = eval_B(Tv,Tv_ref)
global g
  B = g*(Tv-Tv_ref)./Tv_ref;
end