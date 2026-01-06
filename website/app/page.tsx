"use client";

import Link from "next/link";
import Image from "next/image";
import { ArrowRight, Zap, FlaskConical, LineChart, Code, BookOpen } from "lucide-react";
import { motion } from "framer-motion";

export default function Home() {
  const features = [
    {
      icon: Zap,
      title: "Molecular Generation",
      description: "Genetic algorithm-based system that evolves novel fuel molecules from fragment libraries",
    },
    {
      icon: FlaskConical,
      title: "Blend Optimization",
      description: "Deterministic, fast, lab-ready engine that optimizes fuel blends from component libraries",
    },
    {
      icon: LineChart,
      title: "Multi-Objective",
      description: "Optimizes for energy, oxygen balance, volatility, and combustion potential",
    },
    {
      icon: Code,
      title: "Version Agnostic",
      description: "Documentation and tools that work across engine versions",
    },
  ];

  return (
    <div className="flex flex-col w-full overflow-x-hidden">
      {/* Hero Section */}
      <section className="relative min-h-[90vh] flex items-center justify-center overflow-hidden w-full">
        {/* Background gradient with flame effect */}
        <div className="absolute inset-0 bg-gradient-to-b from-ember-dark-red via-background to-background" />
        <div className="absolute inset-0 bg-[radial-gradient(circle_at_50%_50%,oklch(0.65_0.20_35/0.1),transparent_50%)]" />
        
        <div className="container relative z-10 mx-auto px-4 py-20 text-center">
          <motion.h1
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8 }}
            className="mb-6 text-6xl font-bold tracking-tight md:text-8xl"
          >
            <span className="bg-gradient-to-r from-[oklch(0.75_0.22_50)] via-[oklch(0.65_0.20_35)] to-[oklch(0.60_0.18_25)] bg-clip-text text-transparent flame-flicker">
              Ignition Protocol
            </span>
          </motion.h1>
          <motion.p
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.2 }}
            className="mb-4 text-2xl font-semibold text-foreground md:text-3xl"
          >
            An Evolutionary Engine for Molecular Combustion
          </motion.p>
          <motion.p
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.4 }}
            className="mx-auto mb-12 max-w-2xl text-lg text-muted-foreground md:text-xl"
          >
            A computational framework for discovering, filtering, and operationalising synthetic fuel candidates
            through molecular generation and blend optimization.
          </motion.p>
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.6 }}
            className="flex flex-col items-center justify-center gap-4 sm:flex-row"
          >
            <Link
              href="/demo"
              className="group relative overflow-hidden rounded-lg border-2 border-primary bg-primary/10 px-8 py-4 font-semibold text-primary transition-all hover:bg-primary/20 hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.5)]"
            >
              Try the Demo
              <ArrowRight className="ml-2 inline h-5 w-5 transition-transform group-hover:translate-x-1" />
            </Link>
            <Link
              href="/achievements"
              className="rounded-lg border-2 border-border bg-card px-8 py-4 font-semibold text-foreground transition-all hover:border-primary hover:bg-card/80"
            >
              View Achievements
            </Link>
            <Link
              href="/engine"
              className="rounded-lg border-2 border-border bg-card px-8 py-4 font-semibold text-foreground transition-all hover:border-primary hover:bg-card/80"
            >
              View Documentation
            </Link>
          </motion.div>
        </div>
      </section>

      {/* Statistics Section */}
      <section className="container mx-auto px-4 py-12 w-full max-w-full overflow-x-hidden">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="grid gap-6 md:grid-cols-2 lg:grid-cols-4"
        >
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            whileInView={{ opacity: 1, scale: 1 }}
            viewport={{ once: true }}
            transition={{ duration: 0.5, delay: 0.1 }}
            className="rounded-lg border border-primary/50 bg-primary/5 p-6 text-center"
          >
            <p className="mb-2 text-4xl font-bold text-primary">0.1-0.3s</p>
            <p className="text-sm text-muted-foreground">Optimization Speed</p>
          </motion.div>
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            whileInView={{ opacity: 1, scale: 1 }}
            viewport={{ once: true }}
            transition={{ duration: 0.5, delay: 0.2 }}
            className="rounded-lg border border-primary/50 bg-primary/5 p-6 text-center"
          >
            <p className="mb-2 text-4xl font-bold text-primary">2.7M+</p>
            <p className="text-sm text-muted-foreground">Fragments Processed</p>
          </motion.div>
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            whileInView={{ opacity: 1, scale: 1 }}
            viewport={{ once: true }}
            transition={{ duration: 0.5, delay: 0.3 }}
            className="rounded-lg border border-primary/50 bg-primary/5 p-6 text-center"
          >
            <p className="mb-2 text-4xl font-bold text-primary">100%</p>
            <p className="text-sm text-muted-foreground">Chemical Validity</p>
          </motion.div>
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            whileInView={{ opacity: 1, scale: 1 }}
            viewport={{ once: true }}
            transition={{ duration: 0.5, delay: 0.4 }}
            className="rounded-lg border border-primary/50 bg-primary/5 p-6 text-center"
          >
            <p className="mb-2 text-4xl font-bold text-primary">3</p>
            <p className="text-sm text-muted-foreground">Fuel Types Supported</p>
          </motion.div>
        </motion.div>
      </section>

      {/* Features Grid */}
      <section className="container mx-auto px-4 py-20 w-full max-w-full overflow-x-hidden">
        <motion.h2
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="mb-12 text-center text-4xl font-bold"
        >
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Key Features
          </span>
        </motion.h2>
        <div className="grid gap-6 md:grid-cols-2 lg:grid-cols-4">
          {features.map((feature, index) => {
            const Icon = feature.icon;
            return (
              <motion.div
                key={feature.title}
                initial={{ opacity: 0, y: 30 }}
                whileInView={{ opacity: 1, y: 0 }}
                viewport={{ once: true }}
                transition={{ duration: 0.6, delay: index * 0.1 }}
                whileHover={{ scale: 1.05, transition: { duration: 0.2 } }}
                className="group relative overflow-hidden rounded-lg border border-border bg-card p-6 transition-all hover:border-primary hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.3)]"
              >
                <div className="mb-4 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary transition-transform group-hover:scale-110">
                  <Icon className="h-6 w-6" />
                </div>
                <h3 className="mb-2 text-xl font-semibold">{feature.title}</h3>
                <p className="text-muted-foreground">{feature.description}</p>
              </motion.div>
            );
          })}
        </div>
      </section>

      {/* Visual Showcase */}
      <section className="container mx-auto px-4 py-20 w-full max-w-full overflow-x-hidden">
        <motion.div
          initial={{ opacity: 0, y: 30 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.8 }}
          className="max-w-3xl mx-auto"
        >
          <h2 className="mb-6 text-4xl font-bold text-center">
            <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
              Project Overview
            </span>
          </h2>
          <p className="mb-4 text-lg text-muted-foreground text-center">
            Ignition Protocol combines two complementary systems:
          </p>
          <ul className="mb-6 space-y-2 text-muted-foreground">
            <li className="flex items-start gap-2">
              <span className="mt-1 text-primary">•</span>
              <span><strong className="text-foreground">Molecular Generator:</strong> Genetic algorithm-based system that evolves novel fuel molecules</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="mt-1 text-primary">•</span>
              <span><strong className="text-foreground">Blend Optimizer:</strong> Deterministic, fast, lab-ready engine for fuel blend optimization</span>
            </li>
          </ul>
          <div className="text-center">
            <Link
              href="/project"
              className="inline-flex items-center gap-2 text-primary hover:underline"
            >
              Learn More <ArrowRight className="h-4 w-4" />
            </Link>
          </div>
        </motion.div>
      </section>

      {/* Comparison Section */}
      <section className="container mx-auto px-4 py-20 w-full max-w-full overflow-x-hidden">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="mb-12 text-center"
        >
          <h2 className="mb-4 text-4xl font-bold">
            <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
              Why Ignition Protocol?
            </span>
          </h2>
          <p className="text-lg text-muted-foreground">
            A modern approach to fuel design that combines speed, reproducibility, and practical application
          </p>
        </motion.div>
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6, delay: 0.2 }}
          className="mx-auto max-w-4xl"
        >
          <div className="overflow-x-auto rounded-lg border border-border bg-card">
            <table className="w-full">
              <thead>
                <tr className="border-b border-border bg-muted/50">
                  <th className="px-6 py-4 text-left font-semibold">Feature</th>
                  <th className="px-6 py-4 text-center font-semibold">Traditional Methods</th>
                  <th className="px-6 py-4 text-center font-semibold text-primary">Ignition Protocol</th>
                </tr>
              </thead>
              <tbody>
                <tr className="border-b border-border/50">
                  <td className="px-6 py-4 font-medium">Optimization Speed</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Minutes to hours</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">0.1-0.3 seconds</td>
                </tr>
                <tr className="border-b border-border/50 bg-muted/20">
                  <td className="px-6 py-4 font-medium">Reproducibility</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Variable, hard to reproduce</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">100% deterministic</td>
                </tr>
                <tr className="border-b border-border/50">
                  <td className="px-6 py-4 font-medium">Database Scale</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Limited by memory</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">2.7M+ fragments (streaming)</td>
                </tr>
                <tr className="border-b border-border/50 bg-muted/20">
                  <td className="px-6 py-4 font-medium">Chemical Validity</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Manual validation</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">100% automated (RDKit)</td>
                </tr>
                <tr className="border-b border-border/50">
                  <td className="px-6 py-4 font-medium">Output Format</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Raw data, requires processing</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">Lab-ready recipes</td>
                </tr>
                <tr className="border-b border-border/50 bg-muted/20">
                  <td className="px-6 py-4 font-medium">Novel Molecule Integration</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Complex, manual</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">Seamless, automated</td>
                </tr>
                <tr className="border-b border-border/50">
                  <td className="px-6 py-4 font-medium">Home Synthesis Support</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Not available</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">Full feasibility assessment</td>
                </tr>
                <tr>
                  <td className="px-6 py-4 font-medium">Database Scalability</td>
                  <td className="px-6 py-4 text-center text-muted-foreground">Limited to memory</td>
                  <td className="px-6 py-4 text-center font-semibold text-primary">GDB-17 ready (billions)</td>
                </tr>
              </tbody>
            </table>
          </div>
        </motion.div>
      </section>

      {/* CTA Section */}
      <section className="container mx-auto px-4 py-20 w-full max-w-full overflow-x-hidden">
        <motion.div
          initial={{ opacity: 0, scale: 0.95 }}
          whileInView={{ opacity: 1, scale: 1 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="relative overflow-hidden rounded-lg border border-primary/50 bg-gradient-to-r from-ember-dark-red/20 via-background to-ember-dark-red/20 p-12 text-center"
        >
          <div className="relative z-10">
            <motion.div
              animate={{ rotate: [0, 5, -5, 0] }}
              transition={{ duration: 2, repeat: Infinity, repeatDelay: 3 }}
            >
              <BookOpen className="mx-auto mb-4 h-12 w-12 text-primary" />
            </motion.div>
            <h2 className="mb-4 text-3xl font-bold">Full Project Book Available</h2>
            <p className="mb-6 text-muted-foreground">
              The complete Ignition Protocol project book is available at the exhibition stand.
              Visit us to learn more about the methodology, results, and future directions.
            </p>
            <Link
              href="/project"
              className="inline-flex items-center gap-2 rounded-lg border-2 border-primary bg-primary/10 px-6 py-3 font-semibold text-primary transition-all hover:bg-primary/20"
            >
              View Highlights <ArrowRight className="h-4 w-4" />
            </Link>
          </div>
        </motion.div>
      </section>
    </div>
  );
}
