# Exposure time calculator for Mt John telescopes and instruments.

This package will give estimates for the signal to noise ratio possible for given instrument setup and conditions, such as Moon phase.

You can pip install mj_etc as:
```bash
pip install git+https://github.com/CheerfulUser/mt_john_etc.git
```

Once installed the etc can be called as follows. In this example we will calculate the time required to reach a SNR of 5 for a 22nd mag object with the 1.8m using the moaR filter on the MOA camera.

```python
from mj_etc import ETC
etc = ETC('moaR','MOA','1.8m',moon_phase='bright')
etc.time_for_snr(5,mag=22,plot=True)
```
![plot](./figs/test_fig.png)

So we can see that a SNR of 5 should be achieved with a 400s exposure.