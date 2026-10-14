import numpy as np
from des import Env, Resource, lognormal


class Lab:
    def __init__(self, env, rng):
        self.env = env
        self.pcr = Resource(env, capacity=1)
        self.gel = Resource(env, capacity=1)
        self.seq = Resource(env, capacity=1)
        self.metrics = {'completed': 0, 'total_cost': 0.0, 'cycle_times': []}
        self.current_sample = {'steps_completed': 0, 'sample_id': 0, 'cost': 0, 'start_time': 0, 'active': False}

    def run_sample_callback(self):
        sample_id = self.current_sample['sample_id']
        cost_so_far = self.current_sample['cost']
        start_time = self.current_sample['start_time']
        steps_completed = self.current_sample['steps_completed']

        if steps_completed < 1:
            # Step 1: Extraction
            extraction_time = lognormal(self.env.rng, median=0.1, sigma=0.3)
            self.env.schedule(extraction_time, self.step2_callback)
            self.current_sample['steps_completed'] = 1
        elif steps_completed == 1:
            # Step 2: PCR
            if self.pcr.request():
                pcr_time = lognormal(self.env.rng, median=0.1, sigma=0.2)
                self.env.schedule(pcr_time, self.step3_callback)
                self.pcr.release()
                self.current_sample['steps_completed'] = 2
                self.current_sample['cost'] += 20
            else:
                # Wait for PCR
                self.env.schedule(0, self.step2_callback)

    def step2_callback(self, **kw):
        current_sample = self.current_sample
        self.current_sample = {'steps_completed': 1, 'sample_id': current_sample['sample_id'],
                               'cost': current_sample['cost'], 'start_time': current_sample['start_time'],
                               'active': True}

        # Step 2: PCR
        if not self.pcr.request():
            # Wait for PCR
            self.env.schedule(0, self.step2_callback)
            return

        pcr_time = lognormal(self.env.rng, median=0.1, sigma=0.2)
        self.env.schedule(pcr_time, self.step3_callback, cost=current_sample['cost'] + 20)

    def step3_callback(self, **kw):
        current_sample = self.current_sample

        # Step 2: PCR
        self.pcr.release()
        pcr_time = lognormal(self.env.rng, median=0.1, sigma=0.2)
        self.env.schedule(pcr_time, self.step4_callback_1, cost=current_sample['cost'])

        # Save state for next step
        self.current_sample = {'steps_completed': 2, 'sample_id': current_sample['sample_id'],
                               'cost': current_sample['cost'], 'start_time': current_sample['start_time'],
                               'active': True}

    def step4_callback_1(self, **kw):
        # Step 3: Foreshore
        if self.gel.request():
            self.env.schedule(0, self.step4_callback_2)
        else:
            # Wait for gel
            self.env.schedule(0, self.step4_callback_1)

    def step4_callback_2(self, **kw):
        current_sample = self.current_sample

        # Step 3: Foreshore
        if not self.gel.request():
            # Wait for gel
            self.env.schedule(0, self.step4_callback_2)
            return

        self.env.schedule(0, self.step4_callback_3, cost=current_sample['cost'])
        self.gel.release()
        self.current_sample = {'steps_completed': 3, 'sample_id': current_sample['sample_id'],
                               'cost': current_sample['cost'], 'start_time': current_sample['start_time'],
                               'active': True}

    def step4_callback_3(self, **kw):
        current_sample = self.current_sample

        # Step 4: Sequencing
        if self.seq.request():
            self.env.schedule(0, self.step5_callback, cost=current_sample['cost'], start_time=current_sample['start_time'])
        else:
            # Wait for seq
            self.env.schedule(0, self.step4_callback_3)

    def step5_callback(self, **kw):
        current_sample = self.current_sample

        # Step 4: Sequencing
        if not self.seq.request():
            # Wait for seq
            self.env.schedule(0, self.step5_callback)
            return

        seq_time = lognormal(self.env.rng, median=1.0, sigma=0.4)
        self.env.schedule(seq_time, self.step6_callback, cost=current_sample['cost'])
        self.seq.release()

    def step6_callback(self, **kw):
        current_sample = self.current_sample

        # Step 4: Sequencing (complete)
        self.seq.release()
        seq_time = lognormal(self.env.rng, median=1.0, sigma=0.4)

        cost = current_sample['cost'] + 100
        self.env.schedule(seq_time, lambda **x: self.step7_callback(cost=cost,
                                                                   sample_id=current_sample['sample_id'],
                                                                   start_time=current_sample['start_time']))

    def step7_callback(self, **kw):
        cost = kw['cost']
        sample_id = kw['sample_id']
        start_time = kw['start_time']

        cost += 100
        cycle_time = self.env.now - start_time

        self.metrics['completed'] += 1
        self.metrics['total_cost'] += cost
        self.metrics['cycle_times'].append(cycle_time)

        self.current_sample = {'active': False}

    def run_sample(self, sample_id):
        self.current_sample = {
            'steps_completed': 0,
            'sample_id': sample_id,
            'cost': 50,
            'start_time': self.env.now,
            'active': True
        }

        extraction_time = lognormal(self.env.rng, median=0.1, sigma=0.3)
        self.env.schedule(extraction_time, self.run_sample_callback)

        while self.current_sample['active']:
            pass

    def metrics(self):
        cycle_times = self.metrics['cycle_times']
        if not cycle_times:
            return {}
        self.metrics['mean_cycle_time'] = np.mean(cycle_times)
        self.metrics['p95_cycle_time'] = np.percentile(cycle_times, 95)
        return self.metrics


def run(seed: int, days: int = 30) -> dict:
    rng = np.random.default_rng(seed)
    env = Env(rng)
    lab = Lab(env, rng)

    # Poisson arrivals
    lambda_rate = 5
    lambda_inv = 1 / lambda_rate

    def arrival_callback(**kw):
        def callback():
            arrival_time = lognormal(env.rng, median=lambda_inv, sigma=lambda_inv)
            env.schedule(arrival_time, arrival_callback)

            sample_id = lab.metrics['completed']
            lab.run_sample(sample_id)

        env.schedule(0.0, callback)

    def first_arrival_callback(**kw):
        arrival_callback()

    env.schedule(0.0, first_arrival_callback)
    env.run(until=days)

    return lab.metrics()


if __name__ == "__main__":
    metrics = run(seed=0)
    print(metrics)