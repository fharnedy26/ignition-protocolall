"use client";

import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";

interface PropertyComparisonProps {
  properties: Record<string, number>;
  baseline?: Record<string, number>;
}

export function PropertyComparison({ properties, baseline }: PropertyComparisonProps) {
  const data = Object.entries(properties).map(([key, value]) => ({
    property: key,
    optimized: value,
    baseline: baseline?.[key] || 0,
  }));

  return (
    <ResponsiveContainer width="100%" height={300}>
      <BarChart data={data}>
        <CartesianGrid strokeDasharray="3 3" stroke="oklch(0.20 0.05 25)" />
        <XAxis
          dataKey="property"
          stroke="oklch(0.83 0 0)"
          style={{ fontSize: "12px" }}
        />
        <YAxis stroke="oklch(0.83 0 0)" style={{ fontSize: "12px" }} />
        <Tooltip
          contentStyle={{
            backgroundColor: "oklch(0.10 0 0)",
            border: "1px solid oklch(0.65 0.20 35)",
            borderRadius: "0.5rem",
          }}
        />
        <Legend />
        <Bar
          dataKey="optimized"
          fill="oklch(0.65 0.20 35)"
          name="Optimized"
        />
        {baseline && (
          <Bar
            dataKey="baseline"
            fill="oklch(0.50 0 0)"
            name="Baseline"
          />
        )}
      </BarChart>
    </ResponsiveContainer>
  );
}

