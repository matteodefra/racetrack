int originaldmax(int s, int sp, int t) {

	int dist = 0;
	int speed;

	if (s <= sp) {
		if ((t - sp + s) % 2 == 0) {
			speed = s + 1;
			for (int j = 0; j < sp - s; ++j) {
				dist += speed;
				speed += 1;
			}
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			speed -= 2;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
		} else {
			speed = s + 1;
			for (int j = 0; j < sp - s; ++j) {
				dist += speed;
				speed += 1;
			}
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			speed -= 1;
			dist += speed;
			speed -= 1;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
		}
	} else {
		if ((t - s + sp) % 2 == 0) {
			speed = s + 1;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			speed -= 2;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			for (int j = 0; j < s - sp; ++j) {
				dist += speed;
				speed -= 1;
			}
		} else {
			speed = s + 1;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			speed -= 1;
			dist += speed;
			speed -= 1;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			for (int j = 0; j < s - sp; ++j) {
				dist += speed;
				speed -= 1;
			}
		}
	}

	return dist;
}



int originaldmin(int s, int sp, int t) {
	int dist = 0;
	int speed;

	if (s <= sp) {
		if ((t - sp + s) % 2 == 0) {
			speed = s - 1;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			speed += 2;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			for (int j = 0; j < sp - s; ++j) {
				dist += speed;
				speed += 1;
			}
		} else {
			speed = s - 1;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			speed += 1;
			dist += speed;
			speed += 1;
			for (int j = 0; j < (t - sp + s) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
			for (int j = 0; j < sp - s; ++j) {
				dist += speed;
				speed += 1;
			}
		}
	} else {
		if ((t - s + sp) % 2 == 0) {
			speed = s - 1;
			for (int j = 0; j < s - sp; ++j) {
				dist += speed;
				speed -= 1;
			}
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			speed += 2;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
		} else {
			speed = s - 1;
			for (int j = 0; j < s - sp; ++j) {
				dist += speed;
				speed -= 1;
			}
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed -= 1;
			}
			speed += 1;
			dist += speed;
			speed += 1;
			for (int j = 0; j < (t - s + sp) / 2; ++j) {
				dist += speed;
				speed += 1;
			}
		}
	}

	return dist;
}


void get_optimal_interval( int s, int sp, int delta, Interval &interval ) {

	int tinit = abs(s-sp);
	int tstar = abs(s)+abs(sp);

	int dmin, dmax;

	bool found = false;

	for (int t=tinit; t<tstar; t++) {
		dmin = originaldmin(s, sp, t);
		dmax = originaldmax(s, sp, t);

		if (dmin <= delta && dmax >= delta && !found) {
			interval.t1 = t;
			found = true;
		} else if ( ( (!(dmin <= delta) && (dmax>=delta)) || (dmin <= delta && !(dmax>=delta)) ) && found ) {
			interval.t2 = t-1;
			break;
		}
	}

	if (found && interval.t2==-1) {
		interval.t3 = interval.t1;
		interval.t1 = -1;
		interval.t2 = -1;
		return; 
	}
	
	found = false;
	int t = tstar;
	while (!found) {
		dmin = originaldmin(s, sp, t);
		dmax = originaldmax(s, sp, t);

		if (dmin <= delta && dmax >= delta) {
			interval.t3 = t;
			found = true;
		} else {
			t = t+1;
		}
	}

	return;
}