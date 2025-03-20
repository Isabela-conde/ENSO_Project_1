# ENSO_Project_1


ENSO characteristics depend on the background conditions of the tropical Pacific Ocean and the recent ~40-yr trend towards La Niña-like conditions appears to be changing the character of ENSO cycle variability, favoring a shift towards more CP vs EP type ENSO variations. 

Q: are the changing background conditions favouring a different balance of ENSO feedbacks (e.g., zonal advective vs thermocline) that lead to changes in the relative frequency of EP and CP events? 


Suggested data and methods: Heat budget analysis for detrended and undetrended data from various long ocean reanalysis products to reveal the extent to which background conditions have affected the feedbacks and whether those differences are consistent with a shift in the distribution of CP/EP events.

### Data and Methodology-> check if you all agree with this.

**NB#** should have downloaded data from 1976--2024

We compare how the mean state of the two periods (1976--2000, and 2001--2024) affect the two different types of CP and EP ENSO events (I think we should just look at 1 type of event--El Nino??)/ 

We do so by considering components of zonal, meridonal and vertical advection terms, as follows:


$$
    ADV_X &= -u' \frac{\partial \overline{T}}{\partial x} - \overline{u} \frac{\partial T'}{\partial x}\\
    ADV_Y = &= -v' \frac{\partial \overline{T}}{\partial y} - \overline{v} \frac{\partial T'}{\partial y}\\
    ADV_z = &= -w' \frac{\partial \overline{T}}{\partial w} - \overline{w} \frac{\partial T'}{\partial w}\\
$$

In this instance, by separating velocities in terms of anomaly and climatology we can better understand how the mean the difference in mean-state affects the anomalous velocity.

SST anomalies are largely driven by ENSO, therefore when considering ocean velocities anomalies ($u',v',w'$), take this as dominated by the affect of ENSO.

Therefore, zonal advection of mean SST by anomalous zonal currents $-u' \frac{\partial \overline{T}}{\partial t} $, for instance positive zonal anomalies in the equatorial Pacific indicate a stronger easterly current (if this statement is correct-not sure, maybe driven by an increase of trade winds?) then this would increase mean SST in west and decrease it in the east???? No clue tbh. I know the term is 'zonal advection of mean SST by anomalous zonal currents' but unsure how to interpret it. You can interpret the other $ \frac{\partial \overline{T}}{\partial y}$  terms in a similar manner.

The $\overline{w} \frac{\partial T'}{\partial w}$ (or $x/y$) term represents vertical advection of anomalous subsurface temperature by mean upwelling - i.e Ekman pumping??
Using GODAS monthly means from 1980 to 2025, and considering two 25 year periods. We consider the climatology over each period and compute the and anomolus terms for each ENSO CP/EP events, from these terms we take composite for each type of event to see what advective terms contribute most.

Questions I have in the above text: 
- I think we should just look at 1 type of event--El Nino??
- positive zonal anomalies in the equatorial Pacific indicate a stronger easterly current (if this statement is correct-not sure, maybe driven by an increase of trade winds?
  

Calculation of budget:
Will work on this untill someone tells me otherwise! :)

Notes: 

- Please add any thought to this and if you know the answer to any of the questions please let me know or add!
- Also my jupyter notebook has some math in markdown if anybody is intrested.
