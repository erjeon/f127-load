# make_shell.py.retired

임의의 사슬 수로 속 빈 껍질을 만들려고 한 시도이다. 동작하지 않는다.

사슬을 강체로 보고 구면에 한 점씩 배치하는 방식인데, F127 사슬은
긴 축이 8.9 nm 이고 동봉 껍질에서 이웃 사슬 무게중심 간격은 3.37 nm 이다.
사슬 길이가 간격의 2.6 배이므로 실제 미셀에서 사슬들은 서로 얽혀 있다.
구면에 한 점 한 사슬로 놓을 수 있는 물건이 아니다.

50 사슬로 시도했을 때 원자 쌍 55,000 개 이상이 0.9 A 보다 가까웠다.
각 사슬을 자기 축으로 24 번씩 돌려 가며 덜 겹치는 방향을 골라도
줄지 않았다. 반지름을 키우면 겹침은 없어지지만 그것은 미셀이 아니라
빈 껍질이다.

사슬이 구부러져야 자리가 난다. 제대로 하려면 CHARMM-GUI Micelle Builder
를 쓰거나, coarse-grained 로 자가조립시킨 뒤 backmapping 해야 한다.

남겨 둔 것은 같은 시도를 다시 하지 않기 위해서이다.

## load.py

Placed the solute by carving a cavity into an equilibrated micelle and pushing
the polymer out of the way. The pipeline does not use it; `scripts/place_guest.py`
places into the hollow of a shell template instead, and that is the path every
system in the paper was built with.

The two disagreed on how many molecules fit, which is where the wizard's figure
of 54 came from against the 32 the builder could actually place. Its docstring
carried three runnable examples, so anyone reading it would have run something
the rest of the tool does not do.

Kept because the cavity approach is the more general one and may be worth
finishing: it does not need a pre-made hollow template, so it would work for a
micelle of any aggregation number.
