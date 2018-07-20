# **Model Predictive Control (MPC) Project**

[//]: # (Image References)

[image1]: ./md_images/equations.png "Kinematic Equations"
[image2]: ./md_images/equations_latency.png "Equations with Latency"

### The Model

The project uses the kinematic model to control the vehicle’s acceleration, brake and steering. The kinematic model excludes vehicle dynamics such as tire forces, longitudinal and lateral forces, inertia, gravity, air resistance, drag, mass and geometry of vehicle. 

The kinematic model uses the following equations to predict the next (future) vehicle state:

![alt text][image1]

### Timestep Length and Elapsed Duration (N & dt)

The prediction horizon is expressed by T = N * dt where N is the timestep length and dt is the elapsed duration between timesteps. The timestep length, N, determines the number of variables optimized by MPC, to obtain a low cost vector of control inputs and the timestep duration, dt, is the elapsed time between actuations. 

The final N value is 12 and dt value is 0.1 secs, which is the same duration as the actuator latency. This gives a prediction horizon of 1.2 secs, allowing the vehicle to drive round the track smoothly up to 80 mph.


### Polynomial Fitting and MPC Preprocessing

The waypoints are converted to vehicle coordinates, followed by fitting a 3rd order polynomial (`polyfit()`) to the converted waypoints. The polynomial is then evaluated using `polyeval()` where its coefficients are subsequently used to determine the cross-track error (CTE) and orientation error values.


### Model Predictive Control with Latency

The 100 ms actuator latency is incorporated in the computation of the state of vehicle, before `CppAD::ipopt::solve()` uses them to create the reference trajectory.

![alt text][image2]

### Simulation Result

Here's the [link](https://www.dropbox.com/s/xg5seoib6h2fwbu/MPC.mp4?dl=0) to the video result.