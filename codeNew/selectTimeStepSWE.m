function dt = selectTimeStepSWE(dx, dy, avgDepth, gConst, dt, nStep)
etal2 = 5; % TODO
 uul2 = 5; % TODO
 vvl2 = 5; % TODO
hav = avgDepth + etal2;
c = sqrt(gConst*hav);
umax = c + abs(uul2);
vmax = c + abs(vvl2);
dt1 = min(0.5 ./ (umax ./ dx + vmax ./ dy));
if dt1 < dt % TODO evtl andere Fallunterscheidung
	warning(['Time increment reduction necessary in step ' num2str(nStep) ': old value ' num2str(dt) ', new value = ' num2str(dt1) '.']);
end % if
dt = min([dt1 1.05*dt]);
end % function
