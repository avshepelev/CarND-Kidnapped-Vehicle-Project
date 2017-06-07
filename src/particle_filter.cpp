/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <iostream>
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	cout << "ParticleFilter::init" << endl;

	num_particles = 3; //200
	weights.resize(num_particles);

	std::default_random_engine gen;
	// Create normal (Gaussian) distribution for x
	std::normal_distribution<double> dist_x(x, std[0]);
	// Create normal distributions for y and psi
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle par;
		
		par.id = i;
		par.x = dist_x(gen);
		par.y = dist_y(gen);
		par.theta = dist_theta(gen);
		par.weight = 1;

		particles.push_back(par);
	}

	is_initialized = 1;

	// Debugging
	printParticles(particles);
	printWeights();
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	cout << "ParticleFilter::prediction" << endl;

	std::default_random_engine gen;
	
	// Normal distibutions of mean 0 and standard deviation of std_pos
	// Will be added to predicted positions
	std::normal_distribution<double> dist_x(0, std_pos[0]);
	std::normal_distribution<double> dist_y(0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0, std_pos[2]);

	const double theta_dot_dt = yaw_rate * delta_t;

	for (int i = 0; i < num_particles; ++i) {
		double theta_0 = particles[i].theta;

		if (fabs(yaw_rate) > 0.001) {

			particles[i].x += (velocity / yaw_rate) * (sin(theta_0 + theta_dot_dt) - sin(theta_0)) + dist_x(gen);
			particles[i].y += (velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + theta_dot_dt)) + dist_y(gen);
			particles[i].theta += theta_dot_dt + dist_theta(gen);
		}
		else {
			particles[i].x += velocity * delta_t * cos(theta_0) + dist_x(gen);
			particles[i].y += velocity * delta_t * sin(theta_0) + dist_y(gen);
			//particles[i].theta += dist_theta(gen);
			// Theta remains unchanged as yaw rate is zero
		}
		
	}

	// Debugging
	printParticles(particles);
	printWeights();
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// @param predicted - predicted measurement between particle and all map landmarks in sensor range
	// @param observations - actual landmark measurements gathered from lidar
	// Perform nearest neighbor data association and assign each sensor observation a map landmark id associated with it
	// In other words, each observed measurement gets a landmark ID
	cout << "ParticleFilter::dataAssociation" << endl;

	for(int i = 0; i < observations.size(); i++) {
		observations.at(i).id = predicted.at(0).id;
		double min_dist = dist(observations.at(i).x, observations.at(i).y, predicted.at(0).x, predicted.at(0).y);
		//int min_id = predicted.at(1).id;
		observations.at(i).id = 0;

		for(int j = 0; j < predicted.size(); j++) {

			// calculate distance between i_th observation and j_th prediction
			double x1 = observations.at(i).x;
			double y1 = observations.at(i).y;
			double x2 = predicted.at(j).x;
			double y2 = predicted.at(j).y;

			double d = dist(x1, y1, x2, y2);

			if(d < min_dist)
			{
				min_dist = d;
				//min_id = predicted.at(j).id;
				//observations.at(i).id = predicted.at(j).id;
				observations.at(i).id = j;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	//
	// @param sensor_range - range of the sensor
	// @param std_landmark - landmark measurement uncertainties
	// @param observations - landmark measurements
	// @param map_landmarks - 

	//dataAssociation(predicted, *observations)
	cout << "ParticleFilter::updateWeights" << endl;

	std::vector<LandmarkObs> predicted;
	std::vector<LandmarkObs> obs_transformed;
	weights.clear();

	const double denom = (2.0 * M_PI * std_landmark[0] * std_landmark[1]);

	for (int i=0; i<num_particles; ++i) {
		obs_transformed.clear();
		predicted.clear();

		Particle par = particles.at(i);

		// particles are in MAP coordinates
		// observations are in vehicle coordinates
		// transform observations to map coordinates
		
		//cout << "Observations: " << endl;
		//printLandmarks(observations);

		for (int j=0; j<observations.size(); ++j) {
			LandmarkObs obs = observations.at(j);

			double x = par.x + obs.x * cos(par.theta) - obs.y * sin(par.theta);
			double y = par.y + obs.x * sin(par.theta) + obs.y * cos(par.theta);

			//LandmarkObs transformed = {observations.at(j).id, x, y};
			LandmarkObs transformed = {-1, x, y};
			obs_transformed.push_back(transformed);
		}
		//cout << "Transformed observations:" << endl;
		//printLandmarks(obs_transformed);

		// Remove landmarks that are outside sensor range.
		for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
			Map::single_landmark_s landmark = map_landmarks.landmark_list[i];

			if (dist(landmark.x_f, landmark.y_f, par.x, par.y) < sensor_range) {
		        	LandmarkObs l = {landmark.id_i, landmark.x_f, landmark.y_f};
		        	predicted.push_back(l);
			}
		}

		dataAssociation(predicted, obs_transformed);

		double weight = 1.0;
		
		for (int j=0; j<obs_transformed.size(); ++j) {
			double prob = 0;
			LandmarkObs obs = obs_transformed.at(j);

			Map::single_landmark_s closest_landmark = map_landmarks.landmark_list[obs.id - 1];

			double diff_x2 = (obs.x - closest_landmark.x_f) * (obs.x - closest_landmark.x_f);
			double diff_y2 = (obs.y - closest_landmark.y_f) * (obs.y - closest_landmark.y_f);

			double num = exp(-1 * ((diff_x2 / (2.0 * std_landmark[0] * std_landmark[0])) + (diff_y2 / (2.0 * std_landmark[1] * std_landmark[1]))));

			prob = num / denom;
			cout << "prob: " << prob << endl;

			weight *= prob;
		}

		weights.push_back(weight);
		particles.at(i).weight = weight;

	}

	double weights_sum = accumulate(weights.begin(), weights.end(), 0.0);

	for (int i = 0; i < num_particles; ++i) {
		double norm_w = particles[i].weight / weights_sum;
		particles[i].weight = norm_w;
		weights.at(i) = norm_w;
	}

	// Debugging
	printParticles(particles);
	printWeights();
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	cout << "ParticleFilter::resample" << endl;

	cout << "Old particles" << endl;
	printParticles(particles);

	std::vector<Particle> new_particles;
	//new_particles.resize(num_particles);

	std::default_random_engine gen;
	double max_weight = *std::max_element(weights.begin(), weights.end());
	std::uniform_real_distribution<double> dist(0, 2 * max_weight);

	int index = rand() % (num_particles - 1);
	double beta = 0.0;

	for (int i=0; i<num_particles; ++i) {

		beta += dist(gen);
	    
		while (weights[index] < beta){
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
	}

	particles = new_particles;

	// Debugging
	cout << "Num particles: " << num_particles << endl;
	cout << "New particles" << endl;
	printParticles(new_particles);
	
	printWeights();
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

void ParticleFilter::printParticles(std::vector<Particle> pars) {

	cout << "pars.size(): " << pars.size() << endl;

	for (int i = 0; i < pars.size(); ++i) {
		Particle par;
		par = pars.at(i);
		cout << "ID: " << par.id << " X: " << par.x << " Y: " << par.y << " theta: " << par.theta << " weight: " << par.weight << endl;
	}
}

void ParticleFilter::printWeights() {
	for (int j = 0; j < weights.size(); ++j) {
		cout << "Weight " << j << ": " << weights.at(j) << endl;
	}
}

void ParticleFilter::printLandmarks(std::vector<LandmarkObs> landmarks) {
	for (int j = 0; j < landmarks.size(); ++j) {
		LandmarkObs l;
		l = landmarks.at(j);	
		cout << "Landmark " << l.id << ": " << " X: " << l.x << " Y: " << l.y << endl;
	}
}
