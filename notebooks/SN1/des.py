import heapq
import numpy as np


class Env:
    def __init__(self, rng):
        self.now = 0.0
        self.queue = []
        self.rng = rng

    def schedule(self, delay, callback, **kw):
        heapq.heappush(self.queue, (self.now + delay, callback, kw))

    def run(self, until):
        while self.now < until and self.queue:
            time, callback, kw = heapq.heappop(self.queue)
            self.now = time
            callback(**kw)


class Resource:
    def __init__(self, env, capacity):
        self.env = env
        self.capacity = capacity
        self._queue = []

    def request(self):
        if self.capacity > 0:
            self.capacity -= 1
            return True
        else:
            heapq.heappush(self._queue, self.env.now)
            return self._handle()

    def release(self):
        self.capacity += 1
        self._process_queue()

    def _handle(self):
        if self.capacity > 0 and self._queue:
            self.capacity -= 1
            return True
        return False

    def _process_queue(self):
        while self.capacity > 0 and self._queue:
            heapq.heappop(self._queue)
            return True


def lognormal(rng, median, sigma):
    return median * np.exp(sigma * rng.normal())


if __name__ == "__main__":
    rng = np.random.default_rng(42)

    env = Env(rng)
    resource = Resource(env, capacity=1)

    def first_arrival(**kw):
        def callback():
            while resource.request():
                wait_time = lognormal(rng, median=1.0, sigma=0.3)
                env.schedule(wait_time, arrival_done, name="first")
        env.schedule(lognormal(rng, 0.3, 0.2), callback)

    def arrival_done(name, **kw):
        print(f"Finished request {name}")

    env.schedule(lognormal(rng, 0.3, 0.2), first_arrival)
    env.run(until=10.0)