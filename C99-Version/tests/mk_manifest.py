lines = []
for i in range(100):
    P = 2.0 + (i * 0.07) % 10.0
    depth = 0.005 + (i * 0.0003) % 0.008
    lines.append(f'syn:{P:.4f}:{depth:.6f}:10.0:{1000+i}')
open('/tmp/manifest100.txt', 'w').write('\n'.join(lines) + '\n')
print('manifest written')