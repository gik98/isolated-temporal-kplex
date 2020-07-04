delta_aggregation = 30

timestamps = []

try:
    while True:
        time, a, b, ig1, ig2 = input().split("\t")
        a = int(a)
        b = int(b)
        time = int(time)

        aa = min(a, b)
        bb = max(a, b)
        timestamps.append({"time": time, "a": aa, "b": bb})
except EOFError:
    pass

timestamps.sort(key=lambda x: x["time"])

last_timestamp = -1
last_id = -1

edges = {}

for obj in timestamps:
    if obj["time"] != last_timestamp:
        last_id += 1
        last_timestamp = obj["time"]

    if (obj["a"], obj["b"]) not in edges:
        edges[(obj["a"], obj["b"])] = [last_id]
    else:
        edges[(obj["a"], obj["b"])].append(last_id)

for key, val in edges.items():
    times = []
    interval_begin = val[0]
    last_time_end = val[0]
    for time in val:
        if time - delta_aggregation < last_time_end and last_time_end <= time:
            last_time_end = time
        else:
            times.append((interval_begin, last_time_end))
            interval_begin = time
            last_time_end = time

    for time in times:
        print("{} {} {} {}".format(key[0], key[1], time[0], time[1]))
