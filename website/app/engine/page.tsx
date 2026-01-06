"use client";

import Link from "next/link";
import { ArrowRight, BookOpen, Code, Zap, Settings } from "lucide-react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function EnginePage() {
  const sections = [
    {
      icon: Zap,
      title: "Quick Start",
      description: "Get up and running with the engine in minutes",
      href: "/engine/quick-start",
    },
    {
      icon: Code,
      title: "API Reference",
      description: "Complete API documentation with examples",
      href: "/engine/api",
    },
    {
      icon: BookOpen,
      title: "Examples",
      description: "Real-world usage scenarios and code samples",
      href: "/engine/examples",
    },
  ];

  return (
    <div className="container mx-auto px-4 py-12">
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6 }}
        className="mb-12"
      >
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Fuel Blend Optimization Engine
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          A deterministic, fast, lab-ready blend optimization engine that selects and proportions
          real components into optimized fuel blends subject to fuel-type specifications.
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="mb-12 grid gap-6 md:grid-cols-3"
      >
        {sections.map((section, index) => {
          const Icon = section.icon;
          return (
            <Link key={section.href} href={section.href}>
              <motion.div
                initial={{ opacity: 0, y: 30 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.5, delay: 0.3 + index * 0.1 }}
                whileHover={{ scale: 1.05, y: -5 }}
              >
                <Card className="group h-full transition-all hover:border-primary hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.3)]">
                <CardHeader>
                  <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
                    <Icon className="h-6 w-6" />
                  </div>
                  <CardTitle>{section.title}</CardTitle>
                  <CardDescription>{section.description}</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="flex items-center gap-2 text-primary group-hover:gap-4 transition-all">
                    <span className="text-sm font-medium">Learn more</span>
                    <ArrowRight className="h-4 w-4" />
                  </div>
                </CardContent>
              </Card>
              </motion.div>
            </Link>
          );
        })}
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.6 }}
        className="space-y-8"
      >
        <section>
          <h2 className="mb-4 text-3xl font-bold">Key Features</h2>
          <div className="grid gap-4 md:grid-cols-2">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Settings className="h-5 w-5 text-primary" />
                  Deterministic
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Same seed produces identical results, ensuring reproducibility for lab work.
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Zap className="h-5 w-5 text-primary" />
                  Fast Performance
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Optimizes typical component libraries in 0.1-0.3 seconds.
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>Lab-Ready Output</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Volume fractions ready for mixing, with mL and grams per 1000mL recipes.
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardHeader>
                <CardTitle>Fuel-Type Specific</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Supports gasoline, diesel, and jet fuel specifications with appropriate constraints.
                </p>
              </CardContent>
            </Card>
          </div>
        </section>

        <section>
          <h2 className="mb-4 text-3xl font-bold">Performance Metrics</h2>
          <div className="grid gap-4 md:grid-cols-3">
            <Card className="border-primary/50 bg-primary/5">
              <CardHeader>
                <CardTitle className="text-lg">Optimization Speed</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-primary mb-2">0.1-0.3s</p>
                <p className="text-sm text-muted-foreground">
                  For typical component libraries (50-200 components)
                </p>
              </CardContent>
            </Card>
            <Card className="border-primary/50 bg-primary/5">
              <CardHeader>
                <CardTitle className="text-lg">Database Scale</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-primary mb-2">2.7M+</p>
                <p className="text-sm text-muted-foreground">
                  Fragments processed from GDB-11 database
                </p>
              </CardContent>
            </Card>
            <Card className="border-primary/50 bg-primary/5">
              <CardHeader>
                <CardTitle className="text-lg">Chemical Validity</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-primary mb-2">100%</p>
                <p className="text-sm text-muted-foreground">
                  All molecules validated with RDKit
                </p>
              </CardContent>
            </Card>
          </div>
          <Card className="mt-4">
            <CardHeader>
              <CardTitle>Benchmark Profiles</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid gap-3 md:grid-cols-3">
                <div className="rounded-lg border border-border bg-card p-3">
                  <p className="font-semibold text-primary">Ultra-Snappy</p>
                  <p className="text-sm text-muted-foreground">30-40 candidates, 1-2 seconds</p>
                </div>
                <div className="rounded-lg border border-border bg-card p-3">
                  <p className="font-semibold text-primary">Demo</p>
                  <p className="text-sm text-muted-foreground">200-250 candidates, 6-8 seconds</p>
                </div>
                <div className="rounded-lg border border-border bg-card p-3">
                  <p className="font-semibold text-primary">Classroom</p>
                  <p className="text-sm text-muted-foreground">3,000-5,000 candidates, 45-60 seconds</p>
                </div>
              </div>
            </CardContent>
          </Card>
        </section>

        <section>
          <h2 className="mb-4 text-3xl font-bold">Architecture</h2>
          <Card>
            <CardContent className="pt-6">
              <p className="mb-4 text-muted-foreground">
                The engine uses a deterministic greedy-then-refine algorithm that:
              </p>
              <ul className="space-y-2 text-muted-foreground">
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span>Selects components using multi-objective scoring (RON, cetane, LHV, PMI, cost)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span>Enforces fuel-type specific constraints (volatility limits, octane/cetane requirements, density bounds)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span>Respects component caps, novel limits, and delta budgets</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="mt-1 text-primary">•</span>
                  <span>Outputs optimized blends with complete property analysis</span>
                </li>
              </ul>
            </CardContent>
          </Card>
        </section>
      </motion.div>
    </div>
  );
}

