#!/usr/bin/env k8

"use strict";

/*********************************
 * Command-line argument parsing *
 *********************************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

/********************
 * Simpler File I/O *
 ********************/

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

/****************
 * PSMC related *
 ****************/

function psmc_parse(fn, round, scale) {
	let r = { theta0:-1, rho0:-1, C_pi:-1, ri:-1, a:[] };
	let s = 0;
	for (const line of k8_readline(fn)) {
		let m;
		if ((m = /^RD\t(\d+)/.exec(line)) != null) {
			if (parseInt(m[1]) == round)
				s = 1;
		} else if (s == 1 && (m = /^TR\t(\S+)\t(\S+)/.exec(line)) != null) {
			r.theta0 = parseFloat(m[1]);
			r.rho0 = parseFloat(m[2]);
		} else if (s == 1 && (m = /^RI\t(\S+)/.exec(line)) != null) {
			r.ri = parseFloat(m[1]);
		} else if (s == 1 && (m = /^MM\tC_pi: (\S+)/.exec(line)) != null) {
			r.C_pi = parseFloat(m[1]);
		} else if (s == 1 && (m = /^RS\t(\d+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)/.exec(line)) != null) {
			let i = parseInt(m[1]);
			r.a[i] = { t:parseFloat(m[2]), n:parseFloat(m[3]), sigma0:parseFloat(m[4]), sigma1:parseFloat(m[5]), d:0, theta:0 };
		} else if (line == "//") {
			s = 0;
		}
	}
	for (let i = 0; i < r.a.length; ++i) {
		r.a[i].d = r.a[i].t * r.theta0 / scale;
		r.a[i].theta = r.a[i].n * r.theta0 / scale;
	}
	return r;
}

function psmc_select(r, d0) {
	if (d0 <= 0.0) return 0;
	for (let i = 0; i < r.a.length; ++i)
		if (d0 < r.a[i].d)
			return i;
	return -1;
}

function psmc_avg(r, k0, scale) {
	// calculate alpha[]
	let alpha = [], x = 1.0;
	for (let i = k0; i < r.a.length - 1; ++i) {
		alpha.push(x);
		x *= Math.exp(-(r.a[i+1].t - r.a[i].t) / r.a[i].n);
	}
	alpha.push(x);
	alpha.push(0.0);
	// calculate modeled heterozygosity
	let z = 0.0; // C_pi
	for (let i = k0; i < r.a.length; ++i)
		z += r.a[i].n * (alpha[i - k0] - alpha[i - k0 + 1]);
	const C_pi = z;
	z = C_pi * r.theta0 / scale;
	// calculate posterior heterozygosity
	let y = 0.0, f = 0.0;
	for (let i = k0; i < r.a.length; ++i) {
		let t;
		if (i == r.a.length - 1) {
			t = r.a[i].n;
		} else {
			let tao = r.a[i+1].t - r.a[i].t;
			let x = tao / r.a[i].n;
			let tmp = 1.0 / x - 1.0 / (Math.exp(x) - 1);
			if (tmp < 0.0) {
				warn(`Warning: approximation failure at interval ${i}: ${tmp} < 0 (x=${x})`);
				tmp = 0.0;
			} else if (tmp >= 1.0) {
				warn(`Warning: approximation failure at interval ${i}: ${tmp} >= 1 (x=${x})`);
				tmp = 1.0;
			}
			t = tao * tmp;
		}
		f += r.a[i].sigma1;
		y += (1.0 - Math.exp(-r.theta0 * (t + (r.a[i].t - r.a[k0].t)) / scale)) * r.a[i].sigma1;
	}
	return [z, y/f];
}

function main(args) {
	let opt = { round:20, d0:1e-5, scale:100 };
	for (const o of getopt(args, "N:d:s:", [])) {
		if (o.opt == "-N") opt.round = parseInt(o.arg);
		else if (o.opt == "-d") opt.d0 = parseFloat(o.arg);
		else if (o.opt == "-s") opt.scale = parseInt(o.arg);
	}
	if (args.length == 0) {
		print("Usage: split-time.js [options] <in1.psmc> [...]");
		print("Options:");
		print(`  -N INT     pick the INT-th interation [${opt.round}]`);
		print(`  -d FLOAT   sequence divergence cutoff [${opt.d0}]`);
		print(`  -s INT     PSMC scaling used for psmcfa generation [${opt.scale}]`);
		print(`Output:`);
		print(`  1. index of the max time interval that is smaller than {-d}`);
		print(`  2. lower bound of the interval`);
		print(`  3. fraction of genome coalesced <{-d} (modeled)`);
		print(`  4. fraction of genome coalesced <{-d} (posterior)`);
		print(`  5. heterozygosity at t=0 (modeled)`);
		print(`  6. heterozygosity at t=0 (posterior)`);
		print(`  7. heterozygosity at t=Col2 (modeled)`);
		print(`  8. heterozygosity at t=Col2 (posterior)`);
		print(`  9. KL distance between modeled and posterior distributions`);
		return 1;
	}
	for (let j = 0; j < args.length; ++j) {
		warn(`Processing ${args[j]}`);
		let r = psmc_parse(args[j], opt.round, opt.scale);
		let k0 = psmc_select(r, opt.d0);
		if (k0 < 0) {
			warn(`Warning: {-d} is beyond the last time interval for file ${args[j]}. Skipped`);
			continue;
		}
		let sigma0 = 0, sigma1 = 0;
		for (let i = 0; i < k0; ++i)
			sigma0 += r.a[i].sigma0, sigma1 += r.a[i].sigma1;
		const [het0_exp, het0_pos] = psmc_avg(r, 0, opt.scale);
		const [hetk_exp, hetk_pos] = psmc_avg(r, k0, opt.scale);
		print(k0, r.a[k0].d.toFixed(9), sigma0.toFixed(4), sigma1.toFixed(4), het0_exp.toFixed(8), het0_pos.toFixed(8), hetk_exp.toFixed(8), hetk_pos.toFixed(8), r.ri.toFixed(4), args[j]);
	}
}

main(arguments);
