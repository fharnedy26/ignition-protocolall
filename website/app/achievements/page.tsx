"use client";

import { motion } from "framer-motion";
import { CheckCircle2, Zap, FlaskConical, Target, Database, Code, Award } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

export default function AchievementsPage() {
  const achievements = [
    {
      icon: Zap,
      title: "Lightning-Fast Optimization",
      value: "0.1-0.3s",
      description: "Optimizes typical component libraries in under a third of a second",
      highlight: "10-100x faster than traditional methods",
    },
    {
      icon: Database,
      title: "Massive Scale Processing",
      value: "2.7M+",
      description: "Fragments processed from GDB-11 database with memory-efficient streaming",
      highlight: "Handles databases that would crash traditional tools",
    },
    {
      icon: CheckCircle2,
      title: "100% Chemical Validity",
      value: "100%",
      description: "All generated molecules validated with RDKit",
      highlight: "Zero invalid structures in output",
    },
    {
      icon: Target,
      title: "Deterministic Results",
      value: "100%",
      description: "Same seed produces identical results every time",
      highlight: "Perfect reproducibility for scientific work",
    },
    {
      icon: FlaskConical,
      title: "Lab-Ready Outputs",
      value: "Complete",
      description: "Volume fractions, mL, and grams per 1000mL recipes",
      highlight: "Ready for immediate experimental validation",
    },
    {
      icon: Code,
      title: "End-to-End Pipeline",
      value: "Unified",
      description: "From fragment libraries to optimized blend recipes",
      highlight: "Generate → Convert → Optimize → Validate",
    },
    {
      icon: Award,
      title: "Competitive Fitness Scores",
      value: "45.45",
      description: "Generated molecules achieve fitness scores competitive with established fuels",
      highlight: "45.45 vs MTBE 46.03 - near-commercial performance",
    },
  ];

  const scientificRigor = [
    {
      title: "Reproducibility",
      description: "Complete parameter capture in RUN.json files ensures every result can be exactly reproduced",
      details: "All optimization parameters, seeds, and configurations are stored for full traceability",
    },
    {
      title: "Multi-Objective Optimization",
      description: "Balances multiple fuel properties simultaneously (RON, cetane, LHV, PMI, cost)",
      details: "Uses weighted scoring with theme-based optimization (high RON, low PMI, high LHV, low cost)",
    },
    {
      title: "Constraint Enforcement",
      description: "Strict fuel-type specifications (gasoline, diesel, jet) with real-world constraints",
      details: "Enforces volatility limits, octane/cetane requirements, density bounds, and component caps",
    },
    {
      title: "Novel Component Integration",
      description: "Seamlessly integrates novel molecules into optimized blends with configurable caps",
      details: "Supports home-synthesizable molecules with synthesis feasibility assessment",
    },
  ];

  const practicalApplications = [
    "Immediate lab testing with complete mixing recipes",
    "Home synthesis mode: molecules optimized for home synthesis feasibility",
    "Synthesis route prediction with difficulty ratings and yield estimates",
    "Stability testing candidates (oxidation, phase separation, pH)",
    "Multiple fuel types (gasoline, diesel, jet) with appropriate constraints",
    "Delta mode for incremental improvements around existing blends",
    "Waste credit economics for novel components",
    "GPU acceleration support for future scaling",
    "Memory-efficient streaming for billion-molecule databases (GDB-17 ready)",
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
            Key Achievements & Results
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Scientific rigor, performance, and practical applications that set Ignition Protocol apart.
        </p>
      </motion.div>

      {/* Key Metrics */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Performance Highlights</h2>
        <div className="grid gap-6 md:grid-cols-2 lg:grid-cols-3">
          {achievements.map((achievement, index) => {
            const Icon = achievement.icon;
            return (
              <motion.div
                key={achievement.title}
                initial={{ opacity: 0, y: 30 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.5, delay: 0.3 + index * 0.1 }}
                whileHover={{ scale: 1.05, y: -5 }}
              >
                <Card className="h-full border-primary/50 bg-primary/5">
                  <CardHeader>
                    <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
                      <Icon className="h-6 w-6" />
                    </div>
                    <CardTitle>{achievement.title}</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="mb-2 text-3xl font-bold text-primary">{achievement.value}</p>
                    <p className="mb-2 text-sm text-muted-foreground">{achievement.description}</p>
                    <p className="text-xs font-semibold text-accent">{achievement.highlight}</p>
                  </CardContent>
                </Card>
              </motion.div>
            );
          })}
        </div>
      </motion.section>

      {/* Scientific Rigor */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.6 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Scientific Rigor</h2>
        <div className="grid gap-4 md:grid-cols-2">
          {scientificRigor.map((item, index) => (
            <motion.div
              key={item.title}
              initial={{ opacity: 0, x: -20 }}
              animate={{ opacity: 1, x: 0 }}
              transition={{ duration: 0.5, delay: 0.7 + index * 0.1 }}
            >
              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <Award className="h-5 w-5 text-primary" />
                    {item.title}
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="mb-2 text-muted-foreground">{item.description}</p>
                  <p className="text-sm text-muted-foreground/80">{item.details}</p>
                </CardContent>
              </Card>
            </motion.div>
          ))}
        </div>
      </motion.section>

      {/* Practical Applications */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 1.0 }}
      >
        <h2 className="mb-6 text-3xl font-bold">Practical Applications</h2>
        <Card>
          <CardContent className="pt-6">
            <div className="grid gap-3 md:grid-cols-2">
              {practicalApplications.map((app, index) => (
                <motion.div
                  key={index}
                  initial={{ opacity: 0, x: -20 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ duration: 0.4, delay: 1.1 + index * 0.05 }}
                  className="flex items-start gap-3"
                >
                  <CheckCircle2 className="mt-1 h-5 w-5 shrink-0 text-primary" />
                  <span className="text-muted-foreground">{app}</span>
                </motion.div>
              ))}
            </div>
          </CardContent>
        </Card>
      </motion.section>

      {/* Comparison Section */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 1.3 }}
        className="mt-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Why Ignition Protocol?</h2>
        <Card>
          <CardContent className="pt-6">
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead>
                  <tr className="border-b border-border">
                    <th className="pb-3 text-left font-semibold">Feature</th>
                    <th className="pb-3 text-center font-semibold">Traditional Methods</th>
                    <th className="pb-3 text-center font-semibold text-primary">Ignition Protocol</th>
                  </tr>
                </thead>
                <tbody className="space-y-2">
                  <tr className="border-b border-border/50">
                    <td className="py-3">Optimization Speed</td>
                    <td className="py-3 text-center text-muted-foreground">Minutes to hours</td>
                    <td className="py-3 text-center font-semibold text-primary">0.1-0.3 seconds</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Reproducibility</td>
                    <td className="py-3 text-center text-muted-foreground">Variable, hard to reproduce</td>
                    <td className="py-3 text-center font-semibold text-primary">100% deterministic</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Database Scale</td>
                    <td className="py-3 text-center text-muted-foreground">Limited by memory</td>
                    <td className="py-3 text-center font-semibold text-primary">2.7M+ fragments (streaming)</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Chemical Validity</td>
                    <td className="py-3 text-center text-muted-foreground">Manual validation</td>
                    <td className="py-3 text-center font-semibold text-primary">100% automated (RDKit)</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Output Format</td>
                    <td className="py-3 text-center text-muted-foreground">Raw data, requires processing</td>
                    <td className="py-3 text-center font-semibold text-primary">Lab-ready recipes (mL, grams)</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Novel Molecule Integration</td>
                    <td className="py-3 text-center text-muted-foreground">Complex, manual</td>
                    <td className="py-3 text-center font-semibold text-primary">Seamless, automated</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Home Synthesis Support</td>
                    <td className="py-3 text-center text-muted-foreground">Not available</td>
                    <td className="py-3 text-center font-semibold text-primary">Full feasibility assessment</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-3">Experimental Validation</td>
                    <td className="py-3 text-center text-muted-foreground">Separate, disconnected</td>
                    <td className="py-3 text-center font-semibold text-primary">Integrated framework planned</td>
                  </tr>
                  <tr>
                    <td className="py-3">Database Scalability</td>
                    <td className="py-3 text-center text-muted-foreground">Limited to memory</td>
                    <td className="py-3 text-center font-semibold text-primary">GDB-17 ready (billions)</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>
      </motion.section>
    </div>
  );
}

