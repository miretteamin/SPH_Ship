#include "simulation.hpp"
#include <cmath>

using namespace cgp;
using namespace std;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness * (rho - rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
	// To do ...
	//  Fill it with laplacian of W_viscosity
	float const r = norm(p_i - p_j);
	return 45.0f / (3.14159f * pow(h, 6.0f)) * (h - r);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
	// To do ...
	//  Fill it with gradient of W_spiky
	// assert norm(p_i-p_j) > 0 && <= h
	// grad(W_spiky)(p_i - p_j) = -45 / (pi * h^6) * (h - norm(p_i - p_j))^2 * (p_i-p_j) / norm(p_i - p_j)
	float const r = norm(p_i - p_j);
	assert_cgp_no_msg(r <= h);
	return -45.0f / (3.14159f * pow(h, 6.0f)) * pow((h - r), 2) * ((p_i - p_j) / r);
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i - p_j);
	assert_cgp_no_msg(r <= h);
	return 315.0 / (64.0 * 3.14159f * std::pow(h, 9)) * std::pow(h * h - r * r, 3.0f);
}


void update_density(numarray<particle_element>& particles, float h, float m)
{
	// To do: Compute the density value (particles[i].rho) at each particle position
	//  rho_i = \sum_j m W_density(pi,pj)
	int const N = particles.size();

	for (int i = 0; i < N; ++i) {
		float rho = 0.0;
		vec3 p_i = particles[i].p;
		for (int j = 0; j < N; ++j) {
			vec3 p_j = particles[j].p;
			if (norm(p_i - p_j) <= h) {
				rho += m * W_density(p_i, p_j, h);
			}
		}


		particles[i].rho = rho; // to be modified
	}
}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
	for (int i = 0; i < N; ++i) {
		particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
	}

}

// Compute the forces and update the acceleration of the particles
void update_force(float rock_mass, bool boat_alive, bool rock_alive, numarray<particle_element>& boat_particles, numarray<particle_element>& rock_particles, numarray<particle_element>& particles, float h, float m, float nu)
{
	const int N = particles.size();
	const int N_rock = rock_particles.size();
	const int N_boat = boat_particles.size();

	// gravity
	for (int i = 0; i < N; ++i)
		particles[i].f = m * vec3{ 0,-9.81f,0 };

	if (rock_alive) {
		for (int i = 0; i < N_rock; ++i)
			rock_particles[i].f = m * vec3{ 0,-9.81f,0 };
	}
	if (boat_alive) {
		for (int i = 0; i < N_boat; ++i)
			boat_particles[i].f = m * vec3{ 0,-9.81f,0 };
	}

	// Collision between water particles
	for (int i = 0; i < N; ++i) {
		vec3 p_i = particles[i].p;
		vec3 F_pressure_i = { 0.0, 0.0, 0.0 };
		vec3 F_viscosity_i = { 0.0, 0.0, 0.0 };
		float rho_i = particles[i].rho;
		float pr_i = particles[i].pressure;
		vec3 v_i = particles[i].v;
		for (int j = 0; j < N; ++j) {
			vec3 p_j = particles[j].p;
			float const r = norm(p_i - p_j);
			if (i != j && r <= h) {
				float pr_j = particles[j].pressure;
				float rho_j = particles[j].rho;
				vec3 v_j = particles[j].v;
				F_pressure_i += m * (pr_j + pr_i) / (2 * rho_j) * W_gradient_pressure(p_i, p_j, h);
				F_viscosity_i += m * (v_j - v_i) / rho_j * W_laplacian_viscosity(p_i, p_j, h);
			}
		}
		F_pressure_i *= -m / rho_i;
		F_viscosity_i *= m * nu;
		particles[i].f += F_pressure_i + F_viscosity_i;

		// Force from rock particles to water particles
		if (rock_alive) {
			F_pressure_i = { 0.0, 0.0, 0.0 };
			for (int j = 0; j < N_rock; ++j) {
				vec3 p_j = rock_particles[j].p;
				vec3 v_j = rock_particles[j].v;
				float const r = norm(p_i - p_j);
				if (r <= h) {
					float pr_j = rock_particles[j].pressure;
					float rho_j = rock_particles[j].rho;
					particles[i].f += norm(v_j) * (v_j)*rock_mass * 1 / r;
					F_pressure_i += m * (pr_j + pr_i) / (2 * rho_j) * W_gradient_pressure(p_i, p_j, h);
				}
			}
			F_pressure_i *= -m / rho_i;
			// particles[i].f += F_pressure_i;
		}

	}
	// Force from water particles to rock particles
	// if (rock_alive) {
	//     for (int i = 0; i < N_rock; i++) {
	//         vec3 p_i = rock_particles[i].p;
	//         for (int j = 0; j < N; ++j) {
	//             vec3 p_j = particles[j].p;
	//             float const r = norm(p_i-p_j);
	//             if (r <= h) {
	//                 float pr_j = particles[j].pressure;
	//                 rock_particles[i].f += 0.005*0.005/rock_mass * (p_i-p_j) * 1/r * pr_j;
	//             }
	//         }
	//     }
	// }
	if (boat_alive) {
		for (int i = 0; i < N_boat; i++) {
			vec3 p_i = boat_particles[i].p;
			vec3 F_pressure_i = { 0.0, 0.0, 0.0 };
			vec3 F_viscosity_i = { 0.0, 0.0, 0.0 };
			float rho_i = boat_particles[i].rho;
			float pr_i = boat_particles[i].pressure;
			vec3 v_i = boat_particles[i].v;

			for (int j = 0; j < N; ++j) {
				vec3 p_j = particles[j].p;
				float const r = norm(p_i - p_j);
				if (i != j && r <= h) {
					float pr_j = particles[j].pressure;
					float rho_j = particles[j].rho;
					vec3 v_j = particles[j].v;
					F_pressure_i += m * (pr_j + pr_i) / (2 * rho_j) * W_gradient_pressure(p_i, p_j, h);
					F_viscosity_i += m * (v_j - v_i) / rho_j * W_laplacian_viscosity(p_i, p_j, h);
				}
			}
			F_pressure_i *= -m / rho_i;
			F_viscosity_i *= m * nu;
			boat_particles[i].f += F_pressure_i + F_viscosity_i;
		}
	}
}

// numarray<particle_element>& boat_particles
void simulate(float rock_mass, bool boat_alive, bool reset_boat, bool rock_alive, bool reset_rock, bool wave_alive, bool reset_wave, float t, float dt, numarray<particle_element>& particles, numarray<particle_element>& boat_particles, numarray<particle_element>& rock_particles, sph_parameters_structure const& sph_parameters)
{
	float const damping = 0.005f;
	int const N = particles.size();
	float const m = sph_parameters.m;
	int const N_rock = rock_particles.size();
	int const N_boat = boat_particles.size();
	const double PI = 3.141592653589793238463;

	// Add these parameters globally or within your simulation class
	const float waveAmplitude = 0.01f; // Amplitude of the wave
	const float waveLength = 0.3f; // Wavelength
	const float waveSpeed = 2.0f; // Speed at which the wave propagates

	if (wave_alive) {
		float waveNumber = 2 * PI / waveLength;
		float omega = waveNumber * waveSpeed;

		for (int k = 0; k < N; ++k) {
			// is surface particle?
			/*if (particles[k].rho < sph_parameters.rho0 * 0.9f) {
				float phase = waveNumber * particles[k].p.x - omega * t;
				particles[k].p.y += waveAmplitude * sin(phase);
			}*/
			if (particles[k].p.y > -0.5) {
				float phase = waveNumber * particles[k].p.x - omega * t;
				particles[k].p.y += waveAmplitude * sin(phase);
			}
			/*float phase = waveNumber * particles[k].p.x - omega * t;
			particles[k].p.y += waveAmplitude * sin(phase);*/
		}
	}

	// If reset rock button clicked, spawn the rock at a 'random' position at set height. Spawn rock particles around the rock
	if (rock_alive && reset_rock) {
		double angleIncrement = 2 * PI / N_rock;
		float centerX = 0.8 * sin(t);
		float centerY = 1;
		for (int i = 0; i < N_rock; i++) {
			double angle = i * angleIncrement;

			// Calculate coordinates of the point
			double x = centerX + 0.1 * std::sin(angle);
			double y = centerY + 0.1 * std::cos(angle);

			rock_particles[i].p = { x,y,0 };
			rock_particles[i].v = { 0,0,0 };
			rock_particles[i].stopped = false;
		}
	}

	// If reset boat button clicked, spawn the boat at a 'random' position at set height. Spawn boat particles around the boat
	if (boat_alive && reset_boat) {
		for (int i = 0; i < N_boat / 2; i++) {
			boat_particles[i].p = { 0.8 * sin(t) + 0.04 * i, 0.950, 0 };
			boat_particles[i].v = { 0,0,0 };
		}
		for (int i = 0; i < N_boat / 2; i++) {
			boat_particles[N_boat / 2 + i].p = { 0.8 * sin(t) + 0.04 * i, 0.975, 0 };
			boat_particles[N_boat / 2 + i].v = { 0,0,0 };
		}
	}


	update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
	update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
	update_force(rock_mass, boat_alive, rock_alive, boat_particles, rock_particles, particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces    	

	// Collision between rock and floor (before velocity update so that rock will stop)
	float const epsilon = 1e-3f;
	if (rock_alive) {
		for (int k = 0; k < N_rock; ++k)
		{
			vec3& p = rock_particles[k].p;
			vec3& v = rock_particles[k].v;

			// small perturbation to avoid alignment
			if (p.y <= -1) {
				p.y = -1;
				v.y = 0.0f;
				for (int k = 0; k < N_rock; ++k) {
					rock_particles[k].v.y = 0.0f;
					rock_particles[k].stopped = true;
				}
			}
			if (p.x < -1) {
				p.x = -1;
				v.x = 0.0f;
				for (int k = 0; k < N_rock; ++k) {
					rock_particles[k].v.x = 0.0f;
				}
			}
			if (p.x > 1) {
				p.x = 1;
				v.x = 0.0f;
				for (int k = 0; k < N_rock; ++k) {
					rock_particles[k].v.x = 0.0f;
				}
			}

			// if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
			// if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
		}
	}

	// Update water particles position and velocity, limiting speed
	for (int k = 0; k < N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1 - damping) * v + dt * f / m;
		if (norm(v) > 10) {
			v *= 10 / norm(v);
		}
		p = p + dt * v;

	}
	// Update rock particles position and velocity, averaging the velocities of rock particles
	if (rock_alive) {
		vec3 rock_v = { 0,0,0 };
		for (int k = 0; k < N_rock; ++k)
		{
			vec3& v = rock_particles[k].v;
			vec3& f = rock_particles[k].f;

			rock_v += (1 - damping) * v + dt * f / m;
		}
		rock_v /= N_rock;
		for (int k = 0; k < N_rock; ++k)
		{
			if (rock_particles[k].stopped == false) {
				rock_particles[k].v = rock_v;
				rock_particles[k].p += dt * rock_v;
			}

		}
	}

	if (boat_alive) {
		vec3 boat_v = { 0,0,0 };
		for (int k = 0; k < N_boat; ++k)
		{
			vec3& v = boat_particles[k].v;
			vec3& f = boat_particles[k].f;

			boat_v += (1 - damping) * v + dt * f / m;
		}
		boat_v /= N_boat;
		for (int k = 0; k < N_boat; ++k)
		{
			boat_particles[k].v = boat_v;
			boat_particles[k].p += dt * boat_v;
		}
	}



	// Collision with water particles and boundaries

	for (int k = 0; k < N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;

		// small perturbation to avoid alignment
		if (p.y < -1) { p.y = -1 + epsilon * rand_uniform();  v.y *= -0.5f; }
		if (p.x < -1) { p.x = -1 + epsilon * rand_uniform();  v.x *= -0.5f; }
		if (p.x > 1) { p.x = 1 - epsilon * rand_uniform();  v.x *= -0.5f; }
	}

	for (int k = 0; k < N_boat; ++k)
	{
		vec3& p = boat_particles[k].p;
		vec3& v = boat_particles[k].v;

		// small perturbation to avoid alignment
		if (p.y < -1) { p.y = -1 + epsilon * rand_uniform();  v.y *= -0.5f; }
		if (p.x < -1) {
			for (int j = 0; j < N_boat; j++) {
				boat_particles[j].v.x *= -0.5f;
			}
			break;
			// p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;
		}
		if (p.x > 1) {
			for (int j = 0; j < N_boat; j++) {
				boat_particles[j].v.x *= -0.5f;
			}
			break;
			// p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;
		}
	}

}