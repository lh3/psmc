var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

var c, fn_cen = null, ws = 50, step = 100, max = 2e-3;
while ((c = getopt(arguments, "c:w:s:m:")) != null)
	if (c == 'c') fn_cen = getopt.arg;
	else if (c == 'w') ws = parseInt(getopt.arg);
	else if (c == 's') step = parseInt(getopt.arg);
	else if (c == 'm') max = parseFloat(getopt.arg);

ws *= 1000; // in kb

if (arguments.length - getopt.ind < 2) {
	print("Usage: k8 dec2img.js [-c centro.bed] [-w 100] [-s 100] <chr.len> <psmc.dec>");
	exit(1);
}

var f, b = new Bytes();

var len_acc = {}, len_sum = 0;
f = new File(arguments[getopt.ind]);
while (f.readline(b) >= 0) {
	var t = b.toString().split("\t");
	len_acc[t[0]] = len_sum;
	len_sum += parseInt(t[1]);
}
f.close();
var size = Math.pow(2, Math.floor(Math.log(Math.sqrt(len_sum / ws + 1)) / Math.log(2) + 1));

f = new File(arguments[getopt.ind+1]);
var last_chr = null, data = [];
data.length = Math.floor(len_sum / ws) + 1;
for (var i = 0; i < data.length; ++i) data[i] = 0;
while (f.readline(b) >= 0) {
	var m, l = b.toString();
	if ((m = /^DC\s(\S+)\s(\d+)\s(\d+)\s\d+\s(\S+)/.exec(l)) != null) {
		var offset = len_acc[m[1]];
		var st1 = (parseInt(m[2]) - 1) * step;
		var en1 = parseInt(m[3]) * step;
		var stw = Math.floor((st1 + offset) / ws);
		var enw = Math.floor((en1 + offset) / ws);
		var t = parseFloat(m[4]) / step;
		if (stw == enw) { // in one window
			data[stw] += (en1 - st1) / ws * t;
		} else {
			data[stw] += ((stw + 1) * ws - st1 - offset) / ws * t;
			for (var i = stw + 1; i < enw; ++i) data[i] += t;
			data[enw] += (en1 + offset - enw * ws) / ws * t;
		}
	}
}
f.close();

if (fn_cen != null) {
	f = new File(fn_cen);
	while (f.readline(b) >= 0) {
		var t = b.toString().split("\t");
		if (len_acc[t[0]] == null) continue;
		var offset = len_acc[t[0]];
		var stw = Math.floor((parseInt(t[1]) + offset) / ws);
		var enw = Math.floor((parseInt(t[2]) + offset) / ws);
		for (var i = stw; i <= enw; ++i) data[i] = -1;
	}
	f.close();
}

b.destroy();

function rot(n, x, y, rx, ry)
{
	if (ry == 0) {
		if (rx == 1)
			x = n - 1 - x, y = n - 1 - y;
		var t = x; x = y; y = t;
	}
	return [x, y];
}

function d2xy(n, d)
{
	var t = d, x = 0, y = 0;
    for (var s = 1; s < n; s <<= 1) {
        var rx = 1 & (t >> 1);
        var ry = 1 & (t ^ rx);
        var r = rot(s, x, y, rx, ry);
		x = r[0], y = r[1];
        x += s * rx;
        y += s * ry;
        t >>= 2;
    }
	return [x, y];
}

var img = [];
for (var i = 0; i < size; ++i) {
	img[i] = [];
	for (var j = 0; j < size; ++j) img[i][j] = 255;
}

if (max == 0) {
	for (var i = 0; i < data.length; ++i)
		max = max > data[i]? max : data[i];
}

for (var i = 0; i < data.length; ++i) {
	var r = d2xy(size, i);
	if (data[i] > max) data[i] = max;
	if (data[i] >= 0)
		img[r[0]][r[1]] = Math.floor(255 * data[i] / max);
}

print("P2");
print(size + " " + size);
print(255);
for (var i = 0; i < size; ++i) {
	var l = img[i][0];
	for (var j = 1; j < size; ++j)
		l += " " + img[i][j];
	print(l);
}
