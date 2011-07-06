#!/usr/bin/env luajit

-- BEGIN of klib routines
function os.getopt(args, ostr)
	local arg, place = nil, 0;
	return function ()
		if place == 0 then -- update scanning pointer
			place = 1
			if #args == 0 or args[1]:sub(1, 1) ~= '-' then place = 0; return nil end
			if #args[1] >= 2 then
				place = place + 1
				if args[1]:sub(2, 2) == '-' then -- found "--"
					place = 0
					table.remove(args, 1);
					return nil;
				end
			end
		end
		local optopt = args[1]:sub(place, place);
		place = place + 1;
		local oli = ostr:find(optopt);
		if optopt == ':' or oli == nil then -- unknown option
			if optopt == '-' then return nil end
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
			return '?';
		end
		oli = oli + 1;
		if ostr:sub(oli, oli) ~= ':' then -- do not need argument
			arg = nil;
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
		else -- need an argument
			if place <= #args[1] then  -- no white space
				arg = args[1]:sub(place);
			else
				table.remove(args, 1);
				if #args == 0 then -- an option requiring argument is the last one
					place = 0;
					if ostr:sub(1, 1) == ':' then return ':' end
					return '?';
				else arg = args[1] end
			end
			table.remove(args, 1);
			place = 0;
		end
		return optopt, arg;
	end
end

function string:split(sep, n)
	local a, start = {}, 1;
	sep = sep or "%s+";
	repeat
		local b, e = self:find(sep, start);
		if b == nil then
			table.insert(a, self:sub(start));
			break
		end
		a[#a+1] = self:sub(start, b - 1);
		start = e + 1;
		if n and #a == n then
			table.insert(a, self:sub(start));
			break
		end
	until start > #self;
	return a;
end

function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

function math.lgamma(z)
	local x;
	x = 0.1659470187408462e-06     / (z+7);
	x = x + 0.9934937113930748e-05 / (z+6);
	x = x - 0.1385710331296526     / (z+5);
	x = x + 12.50734324009056      / (z+4);
	x = x - 176.6150291498386      / (z+3);
	x = x + 771.3234287757674      / (z+2);
	x = x - 1259.139216722289      / (z+1);
	x = x + 676.5203681218835      / z;
	x = x + 0.9999999999995183;
	return math.log(x) - 5.58106146679532777 - z + (z-0.5) * math.log(z+6.5);
end

function math.lbinom(n, m)
	if m == nil then
		local a = {};
		a[0], a[n] = 0, 0;
		local t = math.lgamma(n+1);
		for m = 1, n-1 do a[m] = t - math.lgamma(m+1) - math.lgamma(n-m+1) end
		return a;
	else return math.lgamma(n+1) - math.lgamma(m+1) - math.lgamma(n-m+1) end
end

function math.bernstein_poly(beta)
	local n = #beta - 1;
	local lbc = math.lbinom(n); -- log binomial coefficients
	return function (t)
		assert(t >= 0 and t <= 1);
		if t == 0 then return beta[1] end
		if t == 1 then return beta[n+1] end
		local sum, logt, logt1 = 0, math.log(t), math.log(1-t);
		for i = 0, n do sum = sum + beta[i+1] * math.exp(lbc[i] + i * logt + (n-i) * logt1) end
		return sum;
	end
end
-- END of klib routines

local min_t, max_t, n_intv, ts_time, alpha, first = 4, 57, 100, false, 0.321, nil
for o, a in os.getopt(arg, 'm:n:sa:M:f:') do
	if o == 'm' then min_t = tonumber(a)
	elseif o == 'M' then max_t = tonumber(a)
	elseif o == 'n' then n_intv = tonumber(a)
	elseif o == 's' then ts_time = true
	elseif o == 'a' then alpha = tonumber(a)
	elseif o == 'f' then first = tonumber(a) <= 0 and nil or tonumber(a)
	end
end

function correct(p, alpha)
	return alpha == nil and p or alpha * (math.exp(p/alpha) - 1)
end

local b, b1 = {{}, {}, {}}, {{}, {}, {}}
if first ~= nil then
	b[1][1], b[2][1], b[3][1] = 1e-7, 1e-7, first * 1e-7
end
local fp = io.xopen(arg[1])
for l in fp:lines() do
	if l:sub(1, 2) == 'CT' then
		local t = l:split('\t')
		if tonumber(t[2]) >= min_t and tonumber(t[2]) <= max_t then
			table.insert(b[1], t[7] / (t[3] + t[5]))
			table.insert(b[2], (t[6] - t[7]) / t[5])
			table.insert(b[3], correct(t[4] / t[3], alpha))
			print('AC', b[1][#b[1]], b[2][#b[2]], b[3][#b[3]])
		end
	end
end
fp:close()

if ts_time then b[1], b[2] = b[2], b[1] end

for i = 1, 3 do
	for j = 1, #b[i] - 1 do
		table.insert(b1[i], (b[i][j+1] - b[i][j]) * (#b[i] - 1))
	end
end

local bp, bp1 = {}, {}
for i = 1, 3 do
	bp[i] = math.bernstein_poly(b[i])
	bp1[i] = math.bernstein_poly(b1[i])
end

for k = 0, n_intv do
	t = k / n_intv
	print('RT', t, bp[1](t), bp[2](t), bp[3](t), bp1[1](t), bp1[2](t), bp1[3](t))
end
